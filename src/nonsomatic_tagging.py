import os
import shlex
import hashlib

from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.vcf import VcfReader, VcfWriter, Position
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


class VcfReader_Database(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 direct_open=False,
                 keep_row_str=False,
                 save_header=False):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.variant_dict = defaultdict(Position)
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.header = ""
        self.save_header = save_header

    def read_vcf(self):
        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        if self.direct_open:
            vcf_fp = open(self.vcf_fn)
            vcf_fo = vcf_fp
        else:
            vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
            vcf_fo = vcf_fp.stdout
        for row in vcf_fo:
            columns = row.strip().split()
            if columns[0][0] == "#":
                if self.save_header:
                    self.header += row
                continue
            chromosome, position = columns[0], columns[1]
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue
            reference, alternate = columns[3], columns[4]
            position = int(position)
            row_str = row if self.keep_row_str else False
            key = (chromosome, position, reference, alternate) if self.ctg_name is None else (position, reference, alternate)

            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              ref_base=reference,
                                              alt_base=alternate,
                                              row_str=row_str
                                              )


def nonsomatic_tag(args):
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

    nonsomatic_filter_variant_dict_id_set_contig = input_variant_dict_id_set_contig.copy()

    total_input = 0
    for k, v in input_variant_dict_id_set_contig.items():
        total_input += len(input_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}...".format(ctg_name))
    print("[INFO] Processing in {}: total input pass calls: {}".format(ctg_name, total_input))

    input_keys_list = list(input_variant_dict_id_set_contig.keys())

    input_inter_pon_variant_dict_id_set_contig_list = []
    pon_vcf_header_info_list = []
    if len(pon_vcf_fn_list) > 0:
        for index, pon_vcf_fn in enumerate(pon_vcf_fn_list):
            input_inter_pon_variant_dict_id_set_contig = defaultdict(set)
            pon_vcf_reader = VcfReader_Database(vcf_fn=str(pon_vcf_fn),
                                                  ctg_name=ctg_name,
                                                  keep_row_str=False,
                                                  save_header=True)
            pon_vcf_reader.read_vcf()
            pon_input_variant_dict = pon_vcf_reader.variant_dict

            pon_variant_dict_id_set_contig = defaultdict(set)

            for k, v in pon_input_variant_dict.items():
                if ctg_name is None:
                    if k[0] not in input_keys_list:
                        continue
                    if str2bool(pon_require_allele_matching_list[index]):
                        for i in range(len(v.alternate_bases)):
                            pon_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[i])
                    else:
                        pon_variant_dict_id_set_contig[k[0]].add(str(k[1]))
                else:
                    if str2bool(pon_require_allele_matching_list[index]):
                        for i in range(len(v.alternate_bases)):
                            pon_variant_dict_id_set_contig[ctg_name].add(str(k[0]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[i])
                    else:
                        pon_variant_dict_id_set_contig[ctg_name].add(str(k[0]))

            for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
                if str2bool(pon_require_allele_matching_list[index]):
                    nonsomatic_filter_variant_dict_id_set_contig[k] = nonsomatic_filter_variant_dict_id_set_contig[k] - \
                                                                      pon_variant_dict_id_set_contig[k]
                    input_inter_pon_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k].intersection(
                        pon_variant_dict_id_set_contig[k])
                else:
                    nonsomatic_filter_variant_dict_id_set_contig[k] = {id for id in nonsomatic_filter_variant_dict_id_set_contig[k] if id.split('\t')[0] not in pon_variant_dict_id_set_contig[k]}
                    input_inter_pon_variant_dict_id_set_contig[k] = input_variant_dict_set_contig[k].intersection(
                        pon_variant_dict_id_set_contig[k])

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