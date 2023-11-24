# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import os
import logging
import shlex

sys.path.insert(0, '/autofs/bal34/lchen/home/ClairST/ClairS-TO')

from time import time
from argparse import ArgumentParser, SUPPRESS
from subprocess import run
from math import log, e
from collections import namedtuple

from shared.vcf import VcfWriter
from shared.utils import str2bool, subprocess_popen
import shared.param as param

import numpy as np

logging.basicConfig(format='%(message)s', level=logging.INFO)

ACGT = 'ACGT'

OutputConfig = namedtuple('OutputConfig', [
    'is_show_reference',
    'quality_score_for_pass',
    'pileup',
    'enable_indel_calling'
])
OutputUtilities = namedtuple('OutputUtilities', [
    'print_debug_message',
    'output',
    'output_header',
    'close_opened_files',
    'gen_output_file'
])


def filtration_value_from(quality_score_for_pass, quality_score, is_reference=False, is_variant=False):
    if is_reference:
        return 'RefCall'
    elif is_variant:
        if quality_score_for_pass is None:
            return 'PASS'
        elif quality_score >= float(quality_score_for_pass):
            return 'PASS'

    return "LowQual"


Phred_Trans = (-10 * log(e, 10))

def quality_score_from(probability, int_format=False, use_phred_qual=True):
    p = float(probability)
    if use_phred_qual:
        tmp = max(Phred_Trans * log(((1.0 - p) + 1e-10) / (p + 1e-10)) + 2.0, 0.0)
    else:
        tmp = max(p, 0.0)

    return int(round(tmp, 4)) if int_format else float(round(tmp, 4))


def argmax(l):
    return max(zip(l, range(len(l))))[1]


def decode_acgt_count(alt_dict, ref_base=None, tumor_coverage=None):
    acgt_count = [0, 0, 0, 0]
    for idx, base in enumerate('ACGT'):
        acgt_count[idx] = alt_dict[base] if base in alt_dict else 0

    if ref_base is not None and tumor_coverage is not None:
        ref_base = ref_base[0].upper()
        if ref_base in 'ACGT':
            #update ref base count
            ref_idx = 'ACGT'.index(ref_base)
            acgt_count[ref_idx] = tumor_coverage - sum(acgt_count)

    AU, CU, GU, TU = acgt_count
    return AU, CU, GU, TU


def output_vcf_from_probability(
        chromosome,
        position,
        reference_base,
        tumor_alt_info,
        probabilities_a,
        probabilities_c,
        probabilities_g,
        probabilities_t,
        probabilities_na,
        probabilities_nc,
        probabilities_ng,
        probabilities_nt,
        likelihood_data_info_list,
        output_config=None,
        vcf_writer=None,
        output_file=None
):
    def decode_alt_info(alt_info):
        alt_info = alt_info.rstrip().split('-')
        read_depth = int(alt_info[0])  # alt_info
        indel_str = alt_info[1] if len(alt_info) > 1 else ''
        seqs = indel_str.split(' ')
        alt_info_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}

        # calcualte the read depth if all positions was deletion or insertion
        if read_depth == 0 and len(alt_info_dict) == 1:
            for k, v in alt_info_dict.items():
                if k[0] == 'D' or k[0] == 'I':
                    read_depth = int(v)

        return alt_info_dict, read_depth

    tumor_alt_info_dict, tumor_read_depth = decode_alt_info(tumor_alt_info)

    alternate_base = reference_base

    probabilities_is_a = probabilities_a[1]
    probabilities_is_c = probabilities_c[1]
    probabilities_is_g = probabilities_g[1]
    probabilities_is_t = probabilities_t[1]
    probabilities_is_na = probabilities_na[1]
    probabilities_is_nc = probabilities_nc[1]
    probabilities_is_ng = probabilities_ng[1]
    probabilities_is_nt = probabilities_nt[1]

    likelihood_a_na = likelihood_data_info_list[0]
    likelihood_c_nc = likelihood_data_info_list[1]
    likelihood_g_ng = likelihood_data_info_list[2]
    likelihood_t_nt = likelihood_data_info_list[3]
    likelihood_a_points_with_zero_one = likelihood_data_info_list[4]
    likelihood_na_points_with_zero_one = likelihood_data_info_list[5]
    likelihood_c_points_with_zero_one = likelihood_data_info_list[6]
    likelihood_nc_points_with_zero_one = likelihood_data_info_list[7]
    likelihood_g_points_with_zero_one = likelihood_data_info_list[8]
    likelihood_ng_points_with_zero_one = likelihood_data_info_list[9]
    likelihood_t_points_with_zero_one = likelihood_data_info_list[10]
    likelihood_nt_points_with_zero_one = likelihood_data_info_list[11]

    a_index = np.digitize(probabilities_is_a, likelihood_a_points_with_zero_one) - 1
    na_index = np.digitize(probabilities_is_na, likelihood_na_points_with_zero_one) - 1

    c_index = np.digitize(probabilities_is_c, likelihood_c_points_with_zero_one) - 1
    nc_index = np.digitize(probabilities_is_nc, likelihood_nc_points_with_zero_one) - 1

    g_index = np.digitize(probabilities_is_g, likelihood_g_points_with_zero_one) - 1
    ng_index = np.digitize(probabilities_is_ng, likelihood_ng_points_with_zero_one) - 1

    t_index = np.digitize(probabilities_is_t, likelihood_t_points_with_zero_one) - 1
    nt_index = np.digitize(probabilities_is_nt, likelihood_nt_points_with_zero_one) - 1

    likelihood_a_na_weight = likelihood_a_na[a_index][na_index] + sys.float_info.epsilon
    likelihood_c_nc_weight = likelihood_c_nc[c_index][nc_index] + sys.float_info.epsilon
    likelihood_g_ng_weight = likelihood_g_ng[g_index][ng_index] + sys.float_info.epsilon
    likelihood_t_nt_weight = likelihood_t_nt[t_index][nt_index] + sys.float_info.epsilon

    probs_is_a = likelihood_a_na_weight * (probabilities_is_a / probabilities_is_na)
    probs_is_c = likelihood_c_nc_weight * (probabilities_is_c / probabilities_is_nc)
    probs_is_g = likelihood_g_ng_weight * (probabilities_is_g / probabilities_is_ng)
    probs_is_t = likelihood_t_nt_weight * (probabilities_is_t / probabilities_is_nt)

    probs_is_acgt = np.array([probs_is_a, probs_is_c, probs_is_g, probs_is_t])
    normalized_probs_is_acgt = probs_is_acgt / np.sum(probs_is_acgt)

    print(''.join([str(item).ljust(25) for item in
                     [" ", "p(Aff|Neg)", "p(Neg|Aff)", "p(Aff)", "p(Neg)", "Normalized_p(Aff|Neg)"]]),
          file=output_file)
    print(''.join([str(item).ljust(25) for item in
                     ["A", probs_is_a, likelihood_a_na_weight, probabilities_is_a, probabilities_is_na, normalized_probs_is_acgt[0]]]),
          file=output_file)
    print(''.join([str(item).ljust(25) for item in
                     ["C", probs_is_c, likelihood_c_nc_weight, probabilities_is_c, probabilities_is_nc, normalized_probs_is_acgt[1]]]),
          file=output_file)
    print(''.join([str(item).ljust(25) for item in
                     ["G", probs_is_g, likelihood_g_ng_weight, probabilities_is_g, probabilities_is_ng, normalized_probs_is_acgt[2]]]),
          file=output_file)
    print(''.join([str(item).ljust(25) for item in
                     ["T", probs_is_t, likelihood_t_nt_weight, probabilities_is_t, probabilities_is_nt, normalized_probs_is_acgt[3]]]),
          file=output_file)
    print('\n', file=output_file)

    p_is_acgt_index = np.argmax(normalized_probs_is_acgt)

    p_is_acgt_value = max(normalized_probs_is_acgt)

    is_a = (p_is_acgt_index == 0)
    is_c = (p_is_acgt_index == 1)
    is_g = (p_is_acgt_index == 2)
    is_t = (p_is_acgt_index == 3)

    is_variant = (is_a and reference_base != 'A') or (is_c and reference_base != 'C') or \
                 (is_g and reference_base != 'G') or (is_t and reference_base != 'T')

    is_reference = not is_variant

    def rank_variant_alt(tumor_alt_info_dict, tumor_read_depth):
        support_alt_dict = {}
        for tumor_alt, tumor_count in tumor_alt_info_dict.items():
            if tumor_alt[0] == 'R':
                continue
            tumor_af = tumor_count / float(tumor_read_depth)
            if tumor_af > 0:
                support_alt_dict[tumor_alt] = tumor_af
        if len(support_alt_dict) == 0:
            return "", 0, 0
        alt_type_list = sorted(support_alt_dict.items(), key=lambda x: x[1], reverse=True)
        best_match_alt_list = []
        for i in range(len(alt_type_list)):
            best_match_alt_list.append(alt_type_list[i][0])
        tumor_supported_reads_count_list = []
        for j in range(len(best_match_alt_list)):
            tumor_supported_reads_count_list.append(tumor_alt_info_dict[best_match_alt_list[j]])

        return best_match_alt_list, tumor_supported_reads_count_list

    if is_variant:
        if tumor_read_depth <= 0:
            print("low tumor coverage")
            return
        best_match_alt_list, tumor_supported_reads_count_list = rank_variant_alt(
            tumor_alt_info_dict, tumor_read_depth)

        best_match_alt = best_match_alt_list[0]
        tumor_supported_reads_count = tumor_supported_reads_count_list[0]

        alternate_base_list = []

        for i in range(len(best_match_alt_list)):
            if best_match_alt_list[i][0] == 'X':
                alternate_base_list.append(best_match_alt_list[i][1])

        if best_match_alt == "":
            return
        if best_match_alt[0] == 'X':
            alternate_base = best_match_alt[1]
            is_SNP = True
            is_variant_a = (is_a and 'A' in alternate_base_list)
            is_variant_c = (is_c and 'C' in alternate_base_list)
            is_variant_g = (is_g and 'G' in alternate_base_list)
            is_variant_t = (is_t and 'T' in alternate_base_list)
            if is_a or is_c or is_g or is_t:
                if not (is_variant_a or is_variant_c or is_variant_g or is_variant_t):
                    is_variant = False
                    is_reference = True
        elif best_match_alt[0] == 'I':
            alternate_base = best_match_alt[1:]
            is_INS = True
        elif best_match_alt[0] == 'D':
            alternate_base = reference_base
            reference_base += best_match_alt[2:]

    if (not output_config.is_show_reference and is_reference) or (
            not is_reference and reference_base == alternate_base):
        return

    if reference_base is None or alternate_base is None:
        return

    # discard Indel
    if (len(reference_base) > 1 or len(alternate_base) > 1) and not output_config.enable_indel_calling:
        return

    if output_config.enable_indel_calling:
        if len(reference_base) == 1 and len(alternate_base) == 1:
            return

    def decode_alt_info(alt_info_dict, read_depth):
        alt_type_list = [{}, {}, {}]  # SNP I D
        ref_num, snp_num, ins_num, del_num = 0, 0, 0, 0
        for alt_type, count in alt_info_dict.items():
            count = int(count)
            if alt_type[0] == 'X':
                alt_type_list[0][alt_type[1]] = count
                snp_num += count
            elif alt_type[0] == 'I':
                alt_type_list[1][alt_type[1:]] = count
                ins_num += count
            elif alt_type[0] == 'D':
                alt_type_list[2][alt_type[1:]] = count
                del_num += count
            elif alt_type[0] == 'R':
                ref_num = count

        return alt_type_list, ref_num, snp_num, ins_num, del_num

    tumor_alt_type_list, tumor_ref_num, tumor_snp_num, tumor_ins_num, tumor_del_num = decode_alt_info(
        alt_info_dict=tumor_alt_info_dict, read_depth=tumor_read_depth)

    if is_reference:
        tumor_supported_reads_count = tumor_ref_num
        alternate_base = "."

    tumor_allele_frequency = min((tumor_supported_reads_count / tumor_read_depth) if tumor_read_depth != 0 else 0.0,
                                 1.0)

    # genotype string
    if is_reference:
        genotype_string = '0/0'
    elif is_variant:
        genotype_string = "0/1" if tumor_allele_frequency < 1.0 else '1/1'

    if is_reference and is_a:
        probabilities_is_a = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_a)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_reference and is_c:
        probabilities_is_c = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_c)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_reference and is_g:
        probabilities_is_g = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_g)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_reference and is_t:
        probabilities_is_t = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_t)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_reference and not is_a and not is_c and not is_g and not is_t:
        if reference_base == 'A':
            probabilities_is_a = normalized_probs_is_acgt[0]
            quality_score = quality_score_from(probabilities_is_a)
            filtration_value = filtration_value_from(
                quality_score_for_pass=output_config.quality_score_for_pass,
                quality_score=quality_score,
                is_reference=is_reference,
                is_variant=is_variant
            )
        elif reference_base == 'C':
            probabilities_is_c = normalized_probs_is_acgt[1]
            quality_score = quality_score_from(probabilities_is_c)
            filtration_value = filtration_value_from(
                quality_score_for_pass=output_config.quality_score_for_pass,
                quality_score=quality_score,
                is_reference=is_reference,
                is_variant=is_variant
            )
        elif reference_base == 'G':
            probabilities_is_g = normalized_probs_is_acgt[2]
            quality_score = quality_score_from(probabilities_is_g)
            filtration_value = filtration_value_from(
                quality_score_for_pass=output_config.quality_score_for_pass,
                quality_score=quality_score,
                is_reference=is_reference,
                is_variant=is_variant
            )
        elif reference_base == 'T':
            probabilities_is_t = normalized_probs_is_acgt[3]
            quality_score = quality_score_from(probabilities_is_t)
            filtration_value = filtration_value_from(
                quality_score_for_pass=output_config.quality_score_for_pass,
                quality_score=quality_score,
                is_reference=is_reference,
                is_variant=is_variant
            )
    elif is_variant and is_a:
        probabilities_is_a = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_a)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_variant and is_c:
        probabilities_is_c = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_c)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_variant and is_g:
        probabilities_is_g = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_g)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    elif is_variant and is_t:
        probabilities_is_t = p_is_acgt_value
        quality_score = quality_score_from(probabilities_is_t)
        filtration_value = filtration_value_from(
            quality_score_for_pass=output_config.quality_score_for_pass,
            quality_score=quality_score,
            is_reference=is_reference,
            is_variant=is_variant,
        )
    else:
        print('ERROR')

    information_string = "."

    AU, CU, GU, TU = decode_acgt_count(tumor_alt_type_list[0], reference_base, tumor_read_depth)

    add_ad_tag = True
    AD = None
    if add_ad_tag:
        AD = str(tumor_supported_reads_count) if is_reference else str(tumor_ref_num) + ',' + str(
            tumor_supported_reads_count)

    vcf_writer.write_row(CHROM=chromosome,
                         POS=position,
                         REF=reference_base,
                         ALT=alternate_base,
                         QUAL=quality_score,
                         FILTER=filtration_value,
                         INFO=information_string,
                         GT=genotype_string,
                         DP=tumor_read_depth,
                         AF=tumor_allele_frequency,
                         AD=AD,
                         AU=AU,
                         CU=CU,
                         GU=GU,
                         TU=TU)

def call_variants_from_probability(args):
    output_config = OutputConfig(
        is_show_reference=args.show_ref,
        quality_score_for_pass=args.qual,
        pileup=args.pileup,
        enable_indel_calling=args.enable_indel_calling
    )

    platform = args.platform

    call_fn = args.call_fn
    if call_fn != "PIPE":
        call_dir = os.path.dirname(call_fn)
        if not os.path.exists(call_dir):
            output = run("mkdir -p {}".format(call_dir), shell=True)
    vcf_writer = VcfWriter(vcf_fn=args.call_fn,
                           ref_fn=args.ref_fn,
                           ctg_name=args.ctg_name,
                           show_ref_calls=args.show_ref,
                           sample_name=args.sample_name,
                           )

    logging.info("[INFO] Calling tumor-only somatic variants ...")
    variant_call_start_time = time()
    prediction_path = args.predict_fn

    if prediction_path != "PIPE":
        if not os.path.exists(prediction_path):
            print("[ERROR] Prediction path not found!")
            return
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, prediction_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    # likelihood_data_fn = args.likelihood_matrix_data
    # likelihood_data = np.loadtxt(likelihood_data_fn)

    likelihood_data_fn = os.path.join('/autofs/bal34/lchen/home/ClairST-Output/resources', 'likelihood_matrixs',
                                      'ont_r10_guppy', '10_equal_distance_division.txt')
    likelihood_data = np.loadtxt(likelihood_data_fn)

    likelihood_data_info_list = []

    likelihood_a_na = likelihood_data[:10]
    likelihood_c_nc = likelihood_data[10:20]
    likelihood_g_ng = likelihood_data[20:30]
    likelihood_t_nt = likelihood_data[30:40]

    likelihood_a_points = likelihood_data[40:41].flatten()[:-1]
    likelihood_na_points = likelihood_data[41:42].flatten()[:-1]
    likelihood_c_points = likelihood_data[42:43].flatten()[:-1]
    likelihood_nc_points = likelihood_data[43:44].flatten()[:-1]
    likelihood_g_points = likelihood_data[44:45].flatten()[:-1]
    likelihood_ng_points = likelihood_data[45:46].flatten()[:-1]
    likelihood_t_points = likelihood_data[46:47].flatten()[:-1]
    likelihood_nt_points = likelihood_data[47:48].flatten()[:-1]

    likelihood_a_points_with_zero = np.insert(likelihood_a_points, 0, 0)
    likelihood_na_points_with_zero = np.insert(likelihood_na_points, 0, 0)
    likelihood_a_points_with_zero_one = np.insert(likelihood_a_points_with_zero, len(likelihood_a_points_with_zero),
                                                  1)
    likelihood_na_points_with_zero_one = np.insert(likelihood_na_points_with_zero,
                                                   len(likelihood_na_points_with_zero), 1)

    likelihood_c_points_with_zero = np.insert(likelihood_c_points, 0, 0)
    likelihood_nc_points_with_zero = np.insert(likelihood_nc_points, 0, 0)
    likelihood_c_points_with_zero_one = np.insert(likelihood_c_points_with_zero, len(likelihood_c_points_with_zero),
                                                  1)
    likelihood_nc_points_with_zero_one = np.insert(likelihood_nc_points_with_zero,
                                                   len(likelihood_nc_points_with_zero), 1)

    likelihood_g_points_with_zero = np.insert(likelihood_g_points, 0, 0)
    likelihood_ng_points_with_zero = np.insert(likelihood_ng_points, 0, 0)
    likelihood_g_points_with_zero_one = np.insert(likelihood_g_points_with_zero, len(likelihood_g_points_with_zero),
                                                  1)
    likelihood_ng_points_with_zero_one = np.insert(likelihood_ng_points_with_zero,
                                                   len(likelihood_ng_points_with_zero), 1)

    likelihood_t_points_with_zero = np.insert(likelihood_t_points, 0, 0)
    likelihood_nt_points_with_zero = np.insert(likelihood_nt_points, 0, 0)
    likelihood_t_points_with_zero_one = np.insert(likelihood_t_points_with_zero, len(likelihood_t_points_with_zero),
                                                  1)
    likelihood_nt_points_with_zero_one = np.insert(likelihood_nt_points_with_zero,
                                                   len(likelihood_nt_points_with_zero), 1)

    likelihood_data_info_list.append(likelihood_a_na)
    likelihood_data_info_list.append(likelihood_c_nc)
    likelihood_data_info_list.append(likelihood_g_ng)
    likelihood_data_info_list.append(likelihood_t_nt)
    likelihood_data_info_list.append(likelihood_a_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_na_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_c_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_nc_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_g_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_ng_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_t_points_with_zero_one)
    likelihood_data_info_list.append(likelihood_nt_points_with_zero_one)

    output_fn = '/autofs/bal13/lchen/home/ClairS-TO_results/ont_r10_guppy/posterior_output/posterior_sample_output.txt'
    output_file = open(output_fn, 'w')

    for row_id, row in enumerate(fo):
        row = row.rstrip().split('\t')
        chromosome, position, reference_base, tumor_alt_info, prediction_a, prediction_c, prediction_g, prediction_t, \
            prediction_na, prediction_nc, prediction_ng, prediction_nt = row[:13]
        probabilities_a = [float(item) for item in prediction_a.split()]
        probabilities_c = [float(item) for item in prediction_c.split()]
        probabilities_g = [float(item) for item in prediction_g.split()]
        probabilities_t = [float(item) for item in prediction_t.split()]
        probabilities_na = [float(item) for item in prediction_na.split()]
        probabilities_nc = [float(item) for item in prediction_nc.split()]
        probabilities_ng = [float(item) for item in prediction_ng.split()]
        probabilities_nt = [float(item) for item in prediction_nt.split()]
        output_vcf_from_probability(
            chromosome,
            position,
            reference_base,
            tumor_alt_info,
            probabilities_a,
            probabilities_c,
            probabilities_g,
            probabilities_t,
            probabilities_na,
            probabilities_nc,
            probabilities_ng,
            probabilities_nt,
            likelihood_data_info_list,
            output_config=output_config,
            vcf_writer=vcf_writer,
            output_file=output_file
        )

    logging.info("[INFO] Total time elapsed: %.2f s" % (time() - variant_call_start_time))

    vcf_writer.close()
    # remove file if on variant in output
    if os.path.exists(args.call_fn):
        vcf_file = open(args.call_fn, 'r').readlines()
        if not len(vcf_file):
            os.remove(args.call_fn)
        for row in vcf_file:
            if row[0] != '#':
                return
        logging.info("[INFO] No vcf output in file {}, remove.".format(args.call_fn))
        os.remove(args.call_fn)


def main():
    parser = ArgumentParser(description="Call variants using trained models and tensors of candidate variants")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required if --gvcf is enabled")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--qual', type=int, default=0,
                        help="If set, variants with >=QUAL will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to samtools, samtools version >= 1.10 is required, default: %(default)s")

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    parser.add_argument('--likelihood_matrix_data', type=str, default=None,
                        help="Likelihood matrix data for Bayes' theorem")

    # options for advanced users
    parser.add_argument('--enable_indel_calling', type=str2bool, default=0,
                        help="EXPERIMENTAL: Call Indel variants, default: disabled")

    # options for debug purpose
    parser.add_argument('--predict_fn', type=str, default="PIPE",
                        help="DEBUG: Output network output probabilities for further analysis")

    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    call_variants_from_probability(args)


if __name__ == "__main__":
    main()