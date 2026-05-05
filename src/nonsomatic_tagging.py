import os
import sys
import io
import hashlib
import subprocess
import threading

from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.vcf import VcfReader
from shared.utils import str2bool, str_none, reference_sequence_from

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

def calculate_file_md5(file_path):
    md5_hash = hashlib.md5()

    with open(file_path, 'rb') as file:
        for chunk in iter(lambda: file.read(4096), b''):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()


def insert_after_line(original_str, target_line, insert_str):
    lines = original_str.split('\n')
    for i, line in enumerate(lines):
        if line == target_line:
            lines.insert(i+1, insert_str.rstrip('\n'))
            break
    return '\n'.join(lines)


def _tagging_summary_line(ctg_name, total_input, tagged_union, remain_pass, pon_fns, pon_hit_totals):
    parts = [
        'ctg={}'.format(ctg_name if ctg_name is not None else '.'),
        'total_input_pass_calls={}'.format(total_input),
        'tagged_non_somatic_union={}'.format(tagged_union),
        'remain_pass_calls={}'.format(remain_pass),
        'num_pon_files={}'.format(len(pon_fns)),
    ]
    for i, (fn, nh) in enumerate(zip(pon_fns, pon_hit_totals), start=1):
        parts.append('PoN{}_hits={}'.format(i, nh))
        parts.append('PoN{}_file={}'.format(i, os.path.basename(str(fn))))
    return ';'.join(parts)


def _pon_has_tabix_index(pon_vcf_fn):
    tbi = pon_vcf_fn + '.tbi' if pon_vcf_fn.endswith('.gz') else pon_vcf_fn + '.tbi'
    return os.path.exists(tbi)


def _build_input_ids_by_contig_pos(input_variant_dict_id_set_contig):
    m = defaultdict(lambda: defaultdict(set))
    for contig, ids in input_variant_dict_id_set_contig.items():
        for id_str in ids:
            pos_str = id_str.split('\t')[0]
            m[contig][pos_str].add(id_str)
    return m


def _parse_pon_line(line, restrict_chroms):
    if line.startswith('#'):
        return None
    line = line.strip()
    if not line:
        return None
    parts = line.split('\t', 5)
    if len(parts) < 5:
        return None
    chromosome, position, reference, alternate = parts[0], parts[1], parts[3], parts[4]
    if restrict_chroms is not None and chromosome not in restrict_chroms:
        return None
    return chromosome, int(position), reference, alternate


def _feed_raw_file_to_gzip_stdin(proc, pon_vcf_fn, md5_hash):
    try:
        with open(pon_vcf_fn, 'rb') as f:
            for chunk in iter(lambda: f.read(262144), b''):
                if md5_hash is not None:
                    md5_hash.update(chunk)
                proc.stdin.write(chunk)
    except BrokenPipeError:
        pass
    finally:
        try:
            proc.stdin.close()
        except Exception:
            pass


def _iter_pon_vcf_records_stream(pon_vcf_fn, ctg_name, input_keys_list, md5_hash=None):
    restrict = frozenset(input_keys_list) if ctg_name is None else None
    if pon_vcf_fn.endswith('.vcf'):
        buf = b''
        with open(pon_vcf_fn, 'rb') as f:
            for chunk in iter(lambda: f.read(262144), b''):
                if md5_hash is not None:
                    md5_hash.update(chunk)
                buf += chunk
                while b'\n' in buf:
                    raw, buf = buf.split(b'\n', 1)
                    line = raw.decode('utf-8', errors='replace')
                    rec = _parse_pon_line(line, restrict)
                    if rec is None:
                        continue
                    if ctg_name is not None and rec[0] != ctg_name:
                        continue
                    yield rec
        if buf:
            line = buf.decode('utf-8', errors='replace')
            rec = _parse_pon_line(line, restrict)
            if rec is not None and (ctg_name is None or rec[0] == ctg_name):
                yield rec
        return

    proc = subprocess.Popen(
        ['gzip', '-dc'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    feeder = threading.Thread(
        target=_feed_raw_file_to_gzip_stdin,
        args=(proc, pon_vcf_fn, md5_hash),
        daemon=True,
    )
    feeder.start()
    try:
        text_out = io.TextIOWrapper(
            proc.stdout, encoding='utf-8', errors='replace', newline='')
        try:
            for line in text_out:
                rec = _parse_pon_line(line, restrict)
                if rec is None:
                    continue
                if ctg_name is not None and rec[0] != ctg_name:
                    continue
                yield rec
        finally:
            text_out.close()
    finally:
        feeder.join()
        proc.wait()


def _apply_pon_streaming(
        pon_vcf_fn,
        ctg_name,
        input_keys_list,
        require_allele,
        input_variant_dict_set_contig,
        input_variant_dict_id_set_contig,
        input_ids_by_contig_pos,
        nonsomatic_filter_variant_dict_id_set_contig,
        skip_pon_md5=False):
    input_inter_pon_variant_dict_id_set_contig = defaultdict(set)

    def apply_one(chromosome, position, reference, alternate):
        alts = alternate.split(',') if ',' in alternate else [alternate]
        if require_allele:
            for alt in alts:
                key = str(position) + '\t' + reference + '\t' + alt
                contig = chromosome if ctg_name is None else ctg_name
                if key in input_variant_dict_id_set_contig[contig]:
                    input_inter_pon_variant_dict_id_set_contig[contig].add(key)
                    nonsomatic_filter_variant_dict_id_set_contig[contig].discard(key)
        else:
            pos_str = str(position)
            contig = chromosome if ctg_name is None else ctg_name
            if pos_str not in input_variant_dict_set_contig[contig]:
                return
            input_inter_pon_variant_dict_id_set_contig[contig].add(pos_str)
            for id_str in input_ids_by_contig_pos[contig][pos_str]:
                nonsomatic_filter_variant_dict_id_set_contig[contig].discard(id_str)

    contigs_to_fetch = [ctg_name] if ctg_name is not None else list(input_keys_list)
    use_tabix = str(pon_vcf_fn).endswith('.gz') and _pon_has_tabix_index(str(pon_vcf_fn))
    used_tabix_ok = False

    def run_gzip_stream():
        md5_obj = None if skip_pon_md5 else hashlib.md5()
        for rec in _iter_pon_vcf_records_stream(
                str(pon_vcf_fn), ctg_name, input_keys_list, md5_hash=md5_obj):
            apply_one(*rec)
        if skip_pon_md5:
            return None
        return md5_obj.hexdigest()

    stream_file_md5 = None

    if use_tabix and len(contigs_to_fetch) == 1 and len(input_keys_list) > 0:
        contig = contigs_to_fetch[0]
        snapshot_filter = None
        try:
            proc = subprocess.Popen(
                ['tabix', '-f', str(pon_vcf_fn), contig],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        except (FileNotFoundError, OSError) as e:
            print("[WARNING] tabix not available or failed for {}: {}. Falling back to full-file stream.".format(
                pon_vcf_fn, str(e)), file=sys.stderr)
            stream_file_md5 = run_gzip_stream()
        else:
            snapshot_filter = {k: set(v) for k, v in nonsomatic_filter_variant_dict_id_set_contig.items()}
            for line in proc.stdout:
                rec = _parse_pon_line(line, None)
                if rec is not None:
                    apply_one(*rec)
            proc.wait()
            stderr_out = proc.stderr.read() if proc.stderr else ''
            if proc.returncode != 0 and stderr_out and 'not a valid contig' not in stderr_out.lower():
                print("[WARNING] tabix failed for {} {}: {}. Falling back to full-file stream.".format(
                    pon_vcf_fn, contig, stderr_out[:200]), file=sys.stderr)
                for k in snapshot_filter:
                    nonsomatic_filter_variant_dict_id_set_contig[k] = set(snapshot_filter[k])
                input_inter_pon_variant_dict_id_set_contig.clear()
                stream_file_md5 = run_gzip_stream()
            else:
                used_tabix_ok = True
    else:
        stream_file_md5 = run_gzip_stream()

    if skip_pon_md5:
        md5_slot = '__skipped__'
    elif used_tabix_ok:
        md5_slot = None
    else:
        md5_slot = stream_file_md5

    return input_inter_pon_variant_dict_id_set_contig, used_tabix_ok, md5_slot


def nonsomatic_tag(args):
    try:
        _nonsomatic_tag_impl(args)
    except MemoryError:
        print("\n[ERROR] Out of memory (OOM) during non-somatic tagging. PoN files may be too large for available RAM.", file=sys.stderr)
        print("[ERROR] Consider: 1) --skip_pon_md5 to skip a full-file read per PoN; 2) Fewer/smaller PoN files; 3) System memory; 4) PoN .tbi for tabix mode.", file=sys.stderr)
        sys.exit(1)
    except OSError as e:
        if 'Cannot allocate memory' in str(e) or 'errno 12' in str(e):
            print("\n[ERROR] Out of memory (OOM) during non-somatic tagging: {}".format(e), file=sys.stderr)
            print("[ERROR] Consider: 1) --skip_pon_md5; 2) Fewer/smaller PoNs; 3) Memory; 4) PoN .tbi index.", file=sys.stderr)
            sys.exit(1)
        raise


def _nonsomatic_tag_impl(args):
    nonsomatic_tag_vcf_header_info = ''
    ctg_name = args.ctg_name
    disable_print_nonsomatic_calls = args.disable_print_nonsomatic_calls
    pileup_vcf_fn = args.pileup_vcf_fn
    panel_of_normals = args.panel_of_normals
    pon_list = list(panel_of_normals.split(',')) if panel_of_normals is not None else []
    panel_of_normals_require_allele_matching = args.panel_of_normals_require_allele_matching
    pon_list_require_allele_matching = list(panel_of_normals_require_allele_matching.split(',')) if panel_of_normals_require_allele_matching is not None else []

    if panel_of_normals is not None:
        if len(pon_list) >= 1:
            pon_vcf_fn_list = pon_list
        else:
            pon_vcf_fn_list = []
    else:
        pon_vcf_fn_list = []

    if panel_of_normals_require_allele_matching is not None:
        if len(pon_list_require_allele_matching) >= 1:
            pon_require_allele_matching_list = pon_list_require_allele_matching
        else:
            pon_require_allele_matching_list = []
    else:
        pon_require_allele_matching_list = []

    pileup_output_vcf_fn = args.output_vcf_fn

    input_vcf_reader = VcfReader(vcf_fn=pileup_vcf_fn,
                                 ctg_name=ctg_name,
                                 show_ref=args.show_ref,
                                 keep_row_str=True,
                                 filter_tag=args.input_filter_tag,
                                 save_header=True)
    input_vcf_reader.read_vcf()
    pileup_variant_dict = input_vcf_reader.variant_dict

    input_variant_dict_set_contig = defaultdict(set)
    input_variant_dict_id_set_contig = defaultdict(set)

    for k, v in pileup_variant_dict.items():
        if ctg_name is None:
            if args.show_ref:
                input_variant_dict_set_contig[k[0]].add(str(k[1]))
                input_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            else:
                if v.filter == "PASS":
                    input_variant_dict_set_contig[k[0]].add(str(k[1]))
                    input_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
        else:
            if args.show_ref:
                input_variant_dict_set_contig[ctg_name].add(str(k))
                input_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            else:
                if v.filter == "PASS":
                    input_variant_dict_set_contig[ctg_name].add(str(k))
                    input_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])

    nonsomatic_filter_variant_dict_id_set_contig = {
        k: set(v) for k, v in input_variant_dict_id_set_contig.items()
    }

    total_input = 0
    for k, v in input_variant_dict_id_set_contig.items():
        total_input += len(input_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}...".format(ctg_name))
    print("[INFO] Processing in {}: total input pass calls: {}".format(ctg_name, total_input))

    input_keys_list = list(input_variant_dict_id_set_contig.keys())

    input_ids_by_contig_pos = _build_input_ids_by_contig_pos(input_variant_dict_id_set_contig)

    input_inter_pon_variant_dict_id_set_contig_list = []
    pon_vcf_header_info_list = []
    if len(pon_vcf_fn_list) > 0:
        for index, pon_vcf_fn in enumerate(pon_vcf_fn_list):
            require_allele = str2bool(pon_require_allele_matching_list[index])

            input_inter_pon_variant_dict_id_set_contig, _, pon_md5_slot = _apply_pon_streaming(
                str(pon_vcf_fn),
                ctg_name,
                input_keys_list,
                require_allele,
                input_variant_dict_set_contig,
                input_variant_dict_id_set_contig,
                input_ids_by_contig_pos,
                nonsomatic_filter_variant_dict_id_set_contig,
                skip_pon_md5=args.skip_pon_md5,
            )

            input_inter_pon_variant_dict_id_set_contig_list.append(
                input_inter_pon_variant_dict_id_set_contig)

            total_filter_by_pon = 0

            for k, v in input_inter_pon_variant_dict_id_set_contig.items():
                total_filter_by_pon += len(input_inter_pon_variant_dict_id_set_contig[k])

            print("[INFO] Processing in {}: tagged by {} PoN: {}".format(ctg_name, str(pon_vcf_fn), total_filter_by_pon))

            if pon_md5_slot == '__skipped__':
                pon_vcf_md5 = 'skipped'
            elif pon_md5_slot is not None:
                pon_vcf_md5 = pon_md5_slot
            else:
                pon_vcf_md5 = calculate_file_md5(str(pon_vcf_fn))
            pon_vcf_header_info = '##INFO=<ID=PoN_{},Number=0,Type=Flag,Description="file={},md5={},allele_matching={},non-somatic variant tagged by panel of normals">'.format(index + 1, str(pon_vcf_fn), pon_vcf_md5, pon_require_allele_matching_list[index])
            pon_vcf_header_info_list.append(pon_vcf_header_info)

            nonsomatic_tag_vcf_header_info += (pon_vcf_header_info + '\n')

    input_inter_all_variant_dict_id_set_contig = defaultdict(set)

    for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
        input_inter_all_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k] - nonsomatic_filter_variant_dict_id_set_contig[k]

    total_filter_by_all = 0

    for k, v in input_inter_all_variant_dict_id_set_contig.items():
        total_filter_by_all += len(input_inter_all_variant_dict_id_set_contig[k])

    total_remain_passes = 0

    for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
        total_remain_passes += len(nonsomatic_filter_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}: tagged by all panel of normals: {}, remained pass calls: {}".format(ctg_name, total_filter_by_all, total_remain_passes))

    pon_hit_totals = []
    for idx in range(len(pon_vcf_fn_list)):
        inter = input_inter_pon_variant_dict_id_set_contig_list[idx]
        pon_hit_totals.append(sum(len(inter[c]) for c in inter))

    summary_line = _tagging_summary_line(
        ctg_name, total_input, total_filter_by_all, total_remain_passes,
        pon_vcf_fn_list, pon_hit_totals)
    print("[INFO] NonSomaticTaggingSummary: {}".format(summary_line))

    contigs_order = major_contigs_order + input_keys_list
    contigs_order_list = sorted(input_keys_list, key=lambda x: contigs_order.index(x))

    original_input_vcf_header = input_vcf_reader.header
    last_filter_line = '##FILTER=<ID=RefCall,Description="Reference call">'
    if nonsomatic_tag_vcf_header_info != '':
        new_nonsomatic_tag_vcf_header = insert_after_line(original_input_vcf_header, last_filter_line, nonsomatic_tag_vcf_header_info)
    else:
        new_nonsomatic_tag_vcf_header = original_input_vcf_header

    with open(pileup_output_vcf_fn, 'w') as output:
        output.write(new_nonsomatic_tag_vcf_header)
        for contig in contigs_order_list:
            all_pos_info = sorted(input_variant_dict_id_set_contig[contig], key=lambda x: int(x.split('\t')[0]))
            for pos_info in all_pos_info:
                pos = int(pos_info.split('\t')[0])
                key = (contig, pos) if ctg_name is None else pos
                ori_row_str = pileup_variant_dict[key].row_str
                if not disable_print_nonsomatic_calls:
                    columns = ori_row_str.split('\t')
                    if pos_info in input_inter_all_variant_dict_id_set_contig[contig]:
                        columns[6] = "NonSomatic"
                        info_str_ori = columns[7]
                        info_str_db = ""
                        if len(input_inter_pon_variant_dict_id_set_contig_list) > 0:
                            for index, input_inter_pon_variant_dict_id_set_contig in enumerate(input_inter_pon_variant_dict_id_set_contig_list):
                                if str2bool(pon_require_allele_matching_list[index]):
                                    if pos_info in input_inter_pon_variant_dict_id_set_contig[contig]:
                                        if info_str_db != "":
                                            info_str_db += ";"
                                        info_str_db += "PoN_{}".format(str(index + 1))
                                else:
                                    if pos_info.split('\t')[0] in input_inter_pon_variant_dict_id_set_contig[contig]:
                                        if info_str_db != "":
                                            info_str_db += ";"
                                        info_str_db += "PoN_{}".format(str(index + 1))

                        columns[7] = info_str_ori + ";" + info_str_db
                    row_str = '\t'.join(columns)
                    output.write(row_str)
                else:
                    if pos_info in nonsomatic_filter_variant_dict_id_set_contig[contig]:
                        output.write(ori_row_str)

    output.close()


def main():
    parser = ArgumentParser(description="Non-somatic tagging for pileup data")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--panel_of_normals', type=str, default=None,
                        help="The path of the panel of normals (PoNs) used for tagging non-somatic variants. Split by ',' for multiple PoNs. Default: if not specified, provided 'gnomad.r2.1.af-ge-0.001.sites.vcf.gz', 'dbsnp.b138.non-somatic.sites.vcf.gz', '1000g-pon.sites.vcf.gz', and 'CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz' PoNs will be included")

    parser.add_argument('--panel_of_normals_require_allele_matching', type=str, default=None,
                        help="Whether to require allele matching when using corresponding PoNs. Split by ',' for multiple PoNs. Default: if not specified, using 'True', 'True', 'False' and 'False' for 'gnomad.r2.1.af-ge-0.001.sites.vcf.gz', 'dbsnp.b138.non-somatic.sites.vcf.gz', '1000g-pon.sites.vcf.gz', and 'CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz' PoNs respectively")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="Output VCF file")

    parser.add_argument('--python', type=str, default="python3",
                        help="Absolute path to the 'python3', default: %(default)s")

    parser.add_argument('--threads', type=int, default=4,
                        help="Max #threads to be used")

    parser.add_argument('--input_filter_tag', type=str_none, default=None,
                        help='Filter variants with tag from the input VCF')

    parser.add_argument('--pypy3', type=str, default="pypy3",
                        help="Absolute path of pypy3, pypy3 >= 3.6 is required")

    parser.add_argument('--parallel', type=str, default="parallel",
                        help="Absolute path of parallel, parallel >= 20191122 is required")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    parser.add_argument("--disable_print_nonsomatic_calls", action='store_true',
                        help="Disable print non-somatic calls. Default: enable non-somatic calls printing")

    parser.add_argument('--skip_pon_md5', action='store_true',
                        help="Skip PoN file MD5 (VCF header uses md5=skipped). Removes one full sequential read per PoN when using tabix, and reduces I/O when large PoNs would otherwise be read twice.")

    global args
    args = parser.parse_args()

    nonsomatic_tag(args)


if __name__ == "__main__":
    main()