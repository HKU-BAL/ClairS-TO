from argparse import ArgumentParser

import sys
import os

from argparse import ArgumentParser
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
from shared.vcf import VcfReader
import subprocess
from src.sort_vcf import compress_index_vcf
from numpy import *
try:
    from scipy.stats import binomtest as binomtest
except ModuleNotFoundError:
     from scipy.stats import binom_test as binomtest

import random
from time import time

seed = int(time())
random.seed(seed)


def tag_germline_variant(args):

    input_vcf_fn = args.input_vcf_fn

    tumor_purity_path = args.tumor_purity_ploidy_output_file

    segment_path = args.tumor_cna_output_file

    if not os.path.exists(tumor_purity_path) or not os.path.exists(segment_path):
        print("[WARNING] Verdict can not obtain final results, not applying verdict tagging!")
        return

    tumor_purity = float(open(tumor_purity_path).read().rstrip().split('\n')[1].split('\t')[1])

    if tumor_purity > 0.6:
        print("[WARNING] Tumor purity estimation {} is higher than 0.6, not applying verdict tagging!".format(tumor_purity))
        return

    input_vcf_reader = VcfReader(
            vcf_fn=input_vcf_fn,
            show_ref=True,
            keep_row_str=True,
            keep_af=True,
            skip_genotype=True,
            save_header=True
    )
    input_vcf_reader.read_vcf()
    input_variant_dict = input_vcf_reader.variant_dict

    cna_file = open(segment_path, 'r')

    seg_chr_list = []
    seg_start_list = []
    seg_end_list = []
    cn_major_list = []
    cn_minor_list = []
    for idx, cna in enumerate(cna_file.readlines()):
        if idx == 0:
            continue
        cna_columns = cna.strip().split('\t')
        seg_chr = str(cna_columns[1].strip('"'))
        seg_start = int(cna_columns[2])
        seg_end = int(cna_columns[3])
        cn_major = int(cna_columns[4])
        cn_minor = int(cna_columns[5])
        seg_chr_list.append(seg_chr)
        seg_start_list.append(seg_start)
        seg_end_list.append(seg_end)
        cn_major_list.append(cn_major)
        cn_minor_list.append(cn_minor)

    cna_file.close()
    ALPHA = 0.01

    for k, v in input_variant_dict.items():
        if v.filter != 'PASS':
            continue
        frequency = v.af
        columns = v.row_str.strip().split('\t')
        depth = int(columns[9].strip().split(':')[2])
        ctg_name = str(columns[0])
        pos = int(columns[1])
        for i in range(len(seg_chr_list)):
            seg_chr = seg_chr_list[i]
            seg_start = seg_start_list[i]
            seg_end = seg_end_list[i]
            cn_major = cn_major_list[i]
            cn_minor = cn_minor_list[i]
            if seg_chr == ctg_name and seg_start <= pos <= seg_end:
                p = tumor_purity
                M = cn_minor
                C = cn_major + cn_minor
                if M == 0:
                    M = C - M
                AF_G1 = (p * M + 1 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                AF_S1 = (p * M + 0 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                P_G1 = binomtest(round(depth * frequency), depth, AF_G1).pvalue
                P_S1 = binomtest(round(depth * frequency), depth, AF_S1).pvalue
                if M != C - M:
                    AF_G2 = (p * (C - M) + 1 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                    P_G2 = binomtest(round(depth * frequency), depth, AF_G2).pvalue
                    if C - M != 0:
                        AF_S2 = (p * (C - M) + 0 * (1 - p)) / (p * C + 2 * (1 - p) + sys.float_info.epsilon)
                        P_S2 = binomtest(round(depth * frequency), depth, AF_S2).pvalue
                    else:
                        AF_S2 = P_S2 = nan
                else:
                    AF_G2 = AF_S2 = P_G2 = P_S2 = nan

                max_prob_germline = max(P_G1, P_G2)
                max_prob_somatic = max(P_S1, P_S2)

                if max_prob_somatic == 0:
                    logodds = inf
                elif max_prob_germline == 0:
                    logodds = -inf
                else:
                    logodds = log10(max_prob_germline) - log10(max_prob_somatic)

                if frequency < 0.05 and 0.2 < p < 0.6:
                    SG_status = 'subclonal somatic'
                    columns[7] += ';Verdict_SubclonalSomatic'

                elif frequency > 0.95:
                    SG_status = 'germline'
                    columns[6] = 'LowQual'
                    columns[7] += ';Verdict_Germline'

                elif max_prob_germline > ALPHA and max_prob_somatic < ALPHA:
                    if logodds < 2:
                        SG_status = 'probable germline'
                    else:
                        if frequency > 0.25:
                            SG_status = 'germline'
                            columns[6] = 'LowQual'
                            columns[7] += ';Verdict_Germline'
                        else:
                            SG_status = 'probable germline'

                elif max_prob_germline < ALPHA and max_prob_somatic > ALPHA:
                    if logodds > -2:
                        SG_status = 'probable somatic'
                    else:
                        SG_status = 'somatic'
                        columns[7] += ';Verdict_Somatic'

                    if nanargmax([P_S1, P_S2]) == 1:
                        M = C - M

                elif max_prob_germline > ALPHA and max_prob_somatic > ALPHA:
                    SG_status = 'ambiguous_both_G_and_S'

                elif max_prob_germline < ALPHA and max_prob_somatic < ALPHA:

                    min_soma_EAF = min(AF_S1, AF_S2)
                    min_germ_EAF = min(AF_G1, AF_G2)

                    if p >= 0.3 and frequency < 0.25 and frequency < min_soma_EAF / 1.5 and min_soma_EAF <= min_germ_EAF:
                        SG_status = 'subclonal somatic'
                        columns[7] += ';Verdict_SubclonalSomatic'

                    elif p >= 0.3 and frequency < 0.25 and frequency < min_germ_EAF / 2.0 and min_germ_EAF < min_soma_EAF:
                        SG_status = 'subclonal somatic'
                        columns[7] += ';Verdict_SubclonalSomatic'

                    elif logodds < -5 and max_prob_somatic > 1e-10:
                        SG_status = 'somatic'
                        columns[7] += ';Verdict_Somatic'

                        if nanargmax([P_S1, P_S2]) == 1:
                            M = C - M

                    elif logodds > 5 and max_prob_germline > 1e-4:
                        SG_status = 'germline'
                        columns[6] = 'LowQual'
                        columns[7] += ';Verdict_Germline'

                    else:
                        SG_status = 'ambiguous_neither_G_nor_S'

                else:
                    SG_status = 'unknown'

                break

        v.row_str = '\t'.join(columns) + '\n'

    with open(args.output_fn, 'w') as f:
        #write header
        header = input_vcf_reader.header
        f.write(header)
        for k, v in input_variant_dict.items():
            f.write(v.row_str)

    if os.path.exists(args.output_fn):
        compress_index_vcf(args.output_fn)


def main():
    parser = ArgumentParser(description="")

    parser.add_argument('--tumor_purity_ploidy_output_file', type=str,
                        default=None,
                        help="Prefix of tumor allele count file")

    parser.add_argument('--tumor_cna_output_file', type=str,
                        default=None,
                        help="Prefix of normal allele count file")

    parser.add_argument('--input_vcf_fn', type=str,
                        default=None,
                        help="Prefix of 1kG alleles file")

    parser.add_argument('--output_fn', type=str,
                        default=None,
                        help="Output file path")

    global args
    args = parser.parse_args()
    tag_germline_variant(args)


if __name__ == "__main__":
    main()