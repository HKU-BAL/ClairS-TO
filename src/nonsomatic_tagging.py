import os
import shlex

from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.vcf import VcfReader, VcfWriter, Position
from shared.utils import str2bool, str_none, reference_sequence_from, subprocess_popen

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


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
            id = str(columns[2])
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue
            reference, alternate = columns[3], columns[4]
            if (len(reference) > 1 or len(alternate) > 1):
                continue
            position = int(position)
            row_str = row if self.keep_row_str else False
            key = (chromosome, position, id) if self.ctg_name is None else (position, id)

            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              ref_base=reference,
                                              alt_base=alternate,
                                              row_str=row_str
                                              )


def nonsomatic_tag(args):

    ctg_name = args.ctg_name
    disable_gnomad_tagging = args.disable_gnomad_tagging
    disable_dbsnp_tagging = args.disable_dbsnp_tagging
    disable_1kgpon_tagging = args.disable_1kgpon_tagging
    disable_print_nonsomatic_calls = args.disable_print_nonsomatic_calls
    pileup_vcf_fn = args.pileup_vcf_fn

    gnomad_vcf_fn = args.gnomad_resource
    dbsnp_vcf_fn = args.dbsnp_resource
    pon_vcf_fn = args.pon_resource
    self_vcf_fn = args.use_own_pon_resource
    pileup_output_vcf_fn = args.output_vcf_fn

    input_vcf_reader = VcfReader(vcf_fn=pileup_vcf_fn,
                                 ctg_name=ctg_name,
                                 show_ref=args.show_ref,
                                 keep_row_str=True,
                                 filter_tag=args.input_filter_tag,
                                 save_header=True)
    input_vcf_reader.read_vcf()
    pileup_variant_dict = input_vcf_reader.variant_dict

    input_variant_dict_id_set_contig = defaultdict(set)

    for k, v in pileup_variant_dict.items():
        if ctg_name is None:
            if args.show_ref:
                input_variant_dict_id_set_contig[k[0]].add(
                    (str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
            else:
                if v.filter == "PASS":
                    input_variant_dict_id_set_contig[k[0]].add(
                        (str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
        else:
            if args.show_ref:
                input_variant_dict_id_set_contig[ctg_name].add(
                    (str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
            else:
                if v.filter == "PASS":
                    input_variant_dict_id_set_contig[ctg_name].add(
                        (str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))

    nonsomatic_filter_variant_dict_id_set_contig = input_variant_dict_id_set_contig.copy()

    input_inter_gnomad_variant_dict_id_set_contig = defaultdict(set)
    input_inter_dbsnp_variant_dict_id_set_contig = defaultdict(set)
    input_inter_pon_variant_dict_id_set_contig = defaultdict(set)
    input_inter_self_variant_dict_id_set_contig = defaultdict(set)

    total_input = 0
    for k, v in input_variant_dict_id_set_contig.items():
        total_input += len(input_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}...".format(ctg_name))
    print("[INFO] Processing in {}: total input pass calls: {}".format(ctg_name, total_input))

    input_keys_list = list(input_variant_dict_id_set_contig.keys())

    if not disable_gnomad_tagging and gnomad_vcf_fn is not None:
        gnomad_vcf_reader = VcfReader(vcf_fn=gnomad_vcf_fn,
                                      ctg_name=ctg_name,
                                      keep_row_str=False,
                                      filter_tag=None,
                                      save_header=True)
        gnomad_vcf_reader.read_vcf()
        gnomad_input_variant_dict = gnomad_vcf_reader.variant_dict

        gnomad_variant_dict_id_set_contig = defaultdict(set)
        for k, v in gnomad_input_variant_dict.items():
            if ctg_name is None:
                if k[0] not in input_keys_list:
                    continue
                gnomad_variant_dict_id_set_contig[k[0]].add(
                    (str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
            else:
                gnomad_variant_dict_id_set_contig[ctg_name].add(
                    (str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))

        for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
            nonsomatic_filter_variant_dict_id_set_contig[k] = nonsomatic_filter_variant_dict_id_set_contig[k] - \
                                                            gnomad_variant_dict_id_set_contig[k]
            input_inter_gnomad_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k].intersection(
                gnomad_variant_dict_id_set_contig[k])

        total_filter_by_gnomad = 0

        for k, v in input_inter_gnomad_variant_dict_id_set_contig.items():
            total_filter_by_gnomad += len(input_inter_gnomad_variant_dict_id_set_contig[k])

        print("[INFO] Processing in {}: tagged by gnomAD resource: {}".format(ctg_name, total_filter_by_gnomad))

    if not disable_dbsnp_tagging and dbsnp_vcf_fn is not None:
        dbsnp_vcf_reader = VcfReader_Database(vcf_fn=dbsnp_vcf_fn,
                                     ctg_name=ctg_name,
                                     keep_row_str=False,
                                     save_header=True)
        dbsnp_vcf_reader.read_vcf()
        dbsnp_input_variant_dict = dbsnp_vcf_reader.variant_dict

        dbsnp_variant_dict_id_set_contig = defaultdict(set)

        for k, v in dbsnp_input_variant_dict.items():
            if ctg_name is None:
                if k[0] not in input_keys_list:
                    continue
                dbsnp_variant_dict_id_set_contig[k[0]].add(
                    (str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
            else:
                dbsnp_variant_dict_id_set_contig[ctg_name].add(
                    (str(k[0]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))

        for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
            nonsomatic_filter_variant_dict_id_set_contig[k] = nonsomatic_filter_variant_dict_id_set_contig[k] - dbsnp_variant_dict_id_set_contig[k]
            input_inter_dbsnp_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k].intersection(dbsnp_variant_dict_id_set_contig[k])

        total_filter_by_dbsnp = 0

        for k, v in input_inter_dbsnp_variant_dict_id_set_contig.items():
            total_filter_by_dbsnp += len(input_inter_dbsnp_variant_dict_id_set_contig[k])

        print("[INFO] Processing in {}: tagged by dbSNP resource: {}".format(ctg_name, total_filter_by_dbsnp))

    if not disable_1kgpon_tagging and pon_vcf_fn is not None:
        pon_vcf_reader = VcfReader(vcf_fn=pon_vcf_fn,
                                   ctg_name=ctg_name,
                                   keep_row_str=False,
                                   filter_tag=None,
                                   save_header=True)
        pon_vcf_reader.read_vcf()
        pon_input_variant_dict = pon_vcf_reader.variant_dict

        pon_variant_dict_id_set_contig = defaultdict(set)

        for k, v in pon_input_variant_dict.items():
            if ctg_name is None:
                pon_variant_dict_id_set_contig[k[0]].add(
                    (str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
            else:
                pon_variant_dict_id_set_contig[ctg_name].add(
                    (str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))

        for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
            nonsomatic_filter_variant_dict_id_set_contig[k] = nonsomatic_filter_variant_dict_id_set_contig[k] - \
                                                            pon_variant_dict_id_set_contig[k]
            input_inter_pon_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k].intersection(
                pon_variant_dict_id_set_contig[k])

        total_filter_by_pon = 0

        for k, v in input_inter_pon_variant_dict_id_set_contig.items():
            total_filter_by_pon += len(input_inter_pon_variant_dict_id_set_contig[k])

        print("[INFO] Processing in {}: tagged by 1000G PoN resource: {}".format(ctg_name, total_filter_by_pon))

    if self_vcf_fn is not None:
        self_vcf_reader = VcfReader(vcf_fn=self_vcf_fn,
                                    ctg_name=ctg_name,
                                    keep_row_str=False,
                                    filter_tag=None,
                                    save_header=True)
        self_vcf_reader.read_vcf()
        self_input_variant_dict = self_vcf_reader.variant_dict

        self_variant_dict_id_set_contig = defaultdict(set)

        for k, v in self_input_variant_dict.items():
            if ctg_name is None:
                self_variant_dict_id_set_contig[k[0]].add(
                    (str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))
            else:
                self_variant_dict_id_set_contig[ctg_name].add(
                    (str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0]))

        for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
            nonsomatic_filter_variant_dict_id_set_contig[k] = nonsomatic_filter_variant_dict_id_set_contig[k] - \
                                                            self_variant_dict_id_set_contig[k]
            input_inter_self_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k].intersection(
                self_variant_dict_id_set_contig[k])

        total_filter_by_self = 0

        for k, v in input_inter_self_variant_dict_id_set_contig.items():
            total_filter_by_self += len(input_inter_self_variant_dict_id_set_contig[k])

        print("[INFO] Processing in  {}: tagged by own PoN resource: {}".format(ctg_name, total_filter_by_self))

    input_inter_all_variant_dict_id_set_contig = defaultdict(set)

    for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
        input_inter_all_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k] - nonsomatic_filter_variant_dict_id_set_contig[k]

    total_filter_by_all = 0

    for k, v in input_inter_all_variant_dict_id_set_contig.items():
        total_filter_by_all += len(input_inter_all_variant_dict_id_set_contig[k])

    total_remain_passes = 0

    for k, v in nonsomatic_filter_variant_dict_id_set_contig.items():
        total_remain_passes += len(nonsomatic_filter_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}: tagged by all database resources: {}, remained pass calls: {}".format(ctg_name, total_filter_by_all, total_remain_passes))

    contigs_order = major_contigs_order + input_keys_list
    contigs_order_list = sorted(input_keys_list, key=lambda x: contigs_order.index(x))

    with open(pileup_output_vcf_fn, 'w') as output:
        output.write(input_vcf_reader.header)
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
                        if pos_info in input_inter_gnomad_variant_dict_id_set_contig[contig]:
                            info_str_db += "gnomAD"
                        if pos_info in input_inter_dbsnp_variant_dict_id_set_contig[contig]:
                            if info_str_db != "":
                                info_str_db += ";"
                            info_str_db += "dbSNP"
                        if pos_info in input_inter_pon_variant_dict_id_set_contig[contig]:
                            if info_str_db != "":
                                info_str_db += ";"
                            info_str_db += "1kGPoN"
                        if pos_info in input_inter_self_variant_dict_id_set_contig[contig]:
                            if info_str_db != "":
                                info_str_db += ";"
                            info_str_db += "ownPoN"
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

    parser.add_argument('--gnomad_resource', type=str, default=None,
                        help="Use gnomAD database resource to tag non-somatic variant")

    parser.add_argument('--dbsnp_resource', type=str, default=None,
                        help="Use dbSNP database resource to tag non-somatic variant")

    parser.add_argument('--pon_resource', type=str, default=None,
                        help="Use 1000G PoN resource to tag non-somatic variant")

    parser.add_argument('--use_own_pon_resource', type=str, default=None,
                        help="Use own variants in VCF format for tagging")

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

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    parser.add_argument(
        "--disable_gnomad_tagging",
        action='store_true',
        help="Disable using gnomAD database for non-somatic variants tagging. Default: enable using gnomAD."
    )

    parser.add_argument(
        "--disable_dbsnp_tagging",
        action='store_true',
        help="Disable using dbSNP database for non-somatic variants tagging. Default: enable using dbSNP."
    )

    parser.add_argument(
        "--disable_1kgpon_tagging",
        action='store_true',
        help="Disable using 1000G PoN for non-somatic variants tagging. Default: enable using 1000G PoN."
    )

    parser.add_argument(
        "--disable_print_nonsomatic_calls",
        action='store_true',
        help="Disable print non-somatic calls. Default: enable non-somatic calls printing."
    )

    global args
    args = parser.parse_args()

    nonsomatic_tag(args)


if __name__ == "__main__":
    main()