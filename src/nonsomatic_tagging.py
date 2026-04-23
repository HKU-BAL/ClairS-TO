import os
import sys
import shlex
import hashlib
import subprocess

from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.vcf import VcfReader
from shared.utils import str2bool, str_none, reference_sequence_from, subprocess_popen

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


def _pon_has_tabix_index(pon_vcf_fn):
    """Check if PoN file has tabix index for region-based streaming."""
    tbi = pon_vcf_fn + '.tbi' if pon_vcf_fn.endswith('.gz') else pon_vcf_fn + '.tbi'
    return os.path.exists(tbi)


def _build_input_ids_by_contig_pos(input_variant_dict_id_set_contig):
    """Input-side index: contig -> position str -> set of full variant id strings (non-allele PoN)."""
    m = defaultdict(lambda: defaultdict(set))
    for contig, ids in input_variant_dict_id_set_contig.items():
        for id_str in ids:
            pos_str = id_str.split('\t')[0]
            m[contig][pos_str].add(id_str)
    return m


def _parse_pon_line(line, restrict_chroms):
    if line.startswith('#'):
        return None
    columns = line.strip().split('\t')
    if len(columns) < 5:
        return None
    chromosome, position, reference, alternate = columns[0], columns[1], columns[3], columns[4]
    if restrict_chroms is not None and chromosome not in restrict_chroms:
        return None
    return chromosome, int(position), reference, alternate


def _iter_pon_vcf_records_gzip(pon_vcf_fn, ctg_name, input_keys_list):
    """Full-file (or uncompressed) PoN stream with same contig filtering as the former full-file PoN load."""
    restrict = frozenset(input_keys_list) if ctg_name is None else None
    if pon_vcf_fn.endswith('.vcf'):
        vcf_fp = open(pon_vcf_fn)
        vcf_fo = vcf_fp
        is_subprocess = False
    else:
        vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (pon_vcf_fn)))
        vcf_fo = vcf_fp.stdout
        is_subprocess = True
    try:
        for line in vcf_fo:
            rec = _parse_pon_line(line, restrict)
            if rec is None:
                continue
            if ctg_name is not None and rec[0] != ctg_name:
                continue
            yield rec
    finally:
        if is_subprocess:
            vcf_fp.stdout.close()
            vcf_fp.wait()


def _apply_pon_streaming(
        pon_vcf_fn,
        ctg_name,
        input_keys_list,
        require_allele,
        input_variant_dict_set_contig,
        input_variant_dict_id_set_contig,
        input_ids_by_contig_pos,
        nonsomatic_filter_variant_dict_id_set_contig):
    """
    Stream PoN once; update nonsomatic_filter and build per-PoN intersection sets.
    Semantics match bulk set intersection / difference on full PoN sets.
    Returns (input_inter_pon_variant_dict_id_set_contig, used_successful_tabix_single_contig).
    """
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
    snapshot_filter = {k: set(v) for k, v in nonsomatic_filter_variant_dict_id_set_contig.items()}
    used_tabix_ok = False

    def run_gzip_stream():
        for rec in _iter_pon_vcf_records_gzip(str(pon_vcf_fn), ctg_name, input_keys_list):
            apply_one(*rec)

    if use_tabix and len(contigs_to_fetch) == 1 and len(input_keys_list) > 0:
        contig = contigs_to_fetch[0]
        try:
            proc = subprocess.Popen(
                ['tabix', '-f', str(pon_vcf_fn), contig],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        except (FileNotFoundError, OSError) as e:
            print("[WARNING] tabix not available or failed for {}: {}. Falling back to full-file stream.".format(
                pon_vcf_fn, str(e)), file=sys.stderr)
            run_gzip_stream()
        else:
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
                run_gzip_stream()
            else:
                used_tabix_ok = True
    else:
        run_gzip_stream()

    return input_inter_pon_variant_dict_id_set_contig, used_tabix_ok


def nonsomatic_tag(args):
    try:
        _nonsomatic_tag_impl(args)
    except MemoryError:
        print("\n[ERROR] Out of memory (OOM) during non-somatic tagging. PoN files may be too large for available RAM.", file=sys.stderr)
        print("[ERROR] Consider: 1) Using fewer/smaller PoN files; 2) Increasing system memory; 3) Ensure PoN files have .tbi index for lower-RAM tabix mode.", file=sys.stderr)
        sys.exit(1)
    except OSError as e:
        if 'Cannot allocate memory' in str(e) or 'errno 12' in str(e):
            print("\n[ERROR] Out of memory (OOM) during non-somatic tagging: {}".format(e), file=sys.stderr)
            print("[ERROR] Consider: 1) Using fewer/smaller PoN files; 2) Increasing system memory; 3) Ensure PoN files have .tbi index.", file=sys.stderr)
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

            input_inter_pon_variant_dict_id_set_contig, used_tabix_ok = _apply_pon_streaming(
                str(pon_vcf_fn),
                ctg_name,
                input_keys_list,
                require_allele,
                input_variant_dict_set_contig,
                input_variant_dict_id_set_contig,
                input_ids_by_contig_pos,
                nonsomatic_filter_variant_dict_id_set_contig,
            )

            if used_tabix_ok:
                print("[INFO] PoN {} loaded by tabix (contig-specific, lower RAM)".format(os.path.basename(pon_vcf_fn)))

            input_inter_pon_variant_dict_id_set_contig_list.append(
                input_inter_pon_variant_dict_id_set_contig)

            total_filter_by_pon = 0

            for k, v in input_inter_pon_variant_dict_id_set_contig.items():
                total_filter_by_pon += len(input_inter_pon_variant_dict_id_set_contig[k])

            print("[INFO] Processing in {}: tagged by {} PoN: {}".format(ctg_name, str(pon_vcf_fn), total_filter_by_pon))

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

    global args
    args = parser.parse_args()

    nonsomatic_tag(args)


if __name__ == "__main__":
    main()