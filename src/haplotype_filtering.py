import os
import shlex
import gc
import subprocess
import concurrent.futures

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import shared.param as param
from shared.vcf import VcfReader, VcfWriter, Position
from shared.utils import str2bool, str_none, reference_sequence_from, subprocess_popen

import sys
import math
from shared.utils import IUPAC_base_to_num_dict as BASE2NUM

LOW_AF_SNV = 0.1
LOW_AF_INDEL = 0.3
min_hom_germline_af = 0.75
eps = 0.5
eps_rse = 0.2
sequence_entropy_threshold = 0.9
flanking = 100

def delete_lines_after(target_str, delimiter):
    lines = target_str.split('\n')
    index = 0
    for i, line in enumerate(lines):
        if delimiter in line:
            index = i
            break
    processed_lines = lines[:index+1]
    processed_str = '\n'.join(processed_lines) + '\n'
    return processed_str


def Binominal(n, k):
    if k > n:
        return 0
    result = 1
    if k > n - k:
        k = n - k
    i = 1
    while i <= k:
        result *= n
        result //= i
        n -= 1
        i += 1

    return result


def fisher_exact(table):
    a, b, c, d = table[0][0], table[0][1], table[1][0], table[1][1]
    if a == b == c == d:
        return 1.0

    p = t = Binominal(a + b, a) * Binominal(c + d, c) / Binominal(a + b + c + d, a + c)

    p_left_side = 0.0
    curP = float(t)
    while (a > 0 and d > 0):
        curP *= a * d
        a -= 1
        b += 1
        c += 1
        d -= 1
        curP /= b * c
        if curP <= t:
            p_left_side += curP

    p += p_left_side

    a, b, c, d = table[0][0], table[0][1], table[1][0], table[1][1]
    p_right_side = 0.0
    curP = float(t)
    while (b > 0 and c > 0):
        curP *= b * c
        a += 1
        b -= 1
        c -= 1
        d += 1
        curP /= a * d
        if curP <=  t:
            p_right_side += curP

    p += p_right_side

    return p


def calculate_sequence_entropy(sequence, entropy_window=None, kmer=5):
    """
    We use a kmer-based sequence entropy calculation to measure the complexity of a region.
    sequence: a chunked sequence around a candidate position, default no_of_positions = flankingBaseNum + 1 + flankingBaseNum
    entropy_window: a maximum entropy window for scanning, if the sequence is larger than the entropy window, a slide
    window would be adopted for measurement.
    kmer: default kmer size for sequence entropy calculation.
    """

    count_repeat_kmer_counts = [0] * (entropy_window + 2)
    count_repeat_kmer_counts[0] = entropy_window

    entropy = [0.0] * (entropy_window + 2)
    for i in range(1, entropy_window + 2):
        e = 1.0 / entropy_window * i
        entropy[i] = e * math.log(e)
    entropy_mul = -1 / math.log(entropy_window)
    entropy_kmer_space = 1 << (2 * kmer)

    kmer_hash_counts = [0] * entropy_kmer_space  # value should smaller than len(seq)
    mask = -1 if kmer > 15 else ~((-1) << (2 * kmer))
    kmer_suffix, kmer_prefix = 0, 0

    i = 0
    i2 = -entropy_window
    entropy_sum = 0.0
    all_entropy_sum = [0.0] * len(sequence)
    while (i2 < len(sequence)):

        if (i < len(sequence)):
            n = BASE2NUM[sequence[i]]
            kmer_suffix = ((kmer_suffix << 2) | n) & mask

            count_repeat_kmer_counts[kmer_hash_counts[kmer_suffix]] -= 1
            entropy_sum -= entropy[kmer_hash_counts[kmer_suffix]]
            kmer_hash_counts[kmer_suffix] += 1
            count_repeat_kmer_counts[kmer_hash_counts[kmer_suffix]] += 1
            entropy_sum += entropy[kmer_hash_counts[kmer_suffix]]

        if i2 >= 0 and i < len(sequence):
            n2 = BASE2NUM[sequence[i2]]
            kmer_prefix = ((kmer_prefix << 2) | n2) & mask  # add base info
            count_repeat_kmer_counts[kmer_hash_counts[kmer_prefix]] -= 1
            entropy_sum -= entropy[kmer_hash_counts[kmer_prefix]]
            kmer_hash_counts[kmer_prefix] -= 1
            count_repeat_kmer_counts[kmer_hash_counts[kmer_prefix]] += 1
            entropy_sum += entropy[kmer_hash_counts[kmer_prefix]]
            all_entropy_sum[i] = entropy_sum
        i += 1
        i2 += 1
    return entropy_sum * entropy_mul


def sqeuence_entropy_from(reference_sequence):
    """
    Calculate sequence entropy in a specific candidate windows, variants in low sequence entropy regions (low
    mappability regions, such as homopolymer, tandem repeat, segmental duplications regions) would more likely have
    more complex variants representation, which is beyond pileup calling. Hence, those candidate variants are re-called by
    full alignment calling.
    We use a kmer-based sequence entropy calculation to measure the complexity of a region, we would directly query the
    chunked reference sequence for sequence entropy calculation for each candidate variant.
    """

    entropy_window = param.no_of_positions
    ref_seq = reference_sequence[flanking - param.flankingBaseNum : flanking + param.flankingBaseNum + 1]
    sequence_entropy = calculate_sequence_entropy(sequence=ref_seq, entropy_window=entropy_window)

    return sequence_entropy


def get_base_list(columns):
    pileup_bases = columns[4]

    base_idx = 0
    base_list = []
    read_end_set = set()
    read_start_set = set()
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]  # add indel seq
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':  # start of read, next base is mq, update mq info
            base_idx += 1
            read_start_set.add(len(base_list) - 1)
        # skip $, the end of read
        if base == "$":
            read_end_set.add(len(base_list) - 1)
        base_idx += 1
    read_start_end_set = read_start_set if len(read_start_set) > len(read_end_set) else read_end_set
    upper_base_counter = Counter([''.join(item).upper() for item in base_list])
    return upper_base_counter, base_list, read_start_end_set


def haplotype_filter_per_pos(args):
    pos = args.pos
    ctg_name = args.ctg_name
    ref_base = args.ref_base
    alt_base = args.alt_base
    tumor_bam_fn = args.tumor_bam_fn
    ref_fn = args.ref_fn
    samtools = args.samtools
    min_bq = args.min_bq
    min_mq = args.min_mq
    debug = args.debug
    max_co_exist_read_num = args.min_alt_coverage
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    is_ins = len(ref_base) == 1 and len(alt_base) > 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    pass_hetero = True
    pass_homo = True
    pass_hetero_both_side = True
    pass_read_start_end = True
    pass_bq = True
    pass_mq = True
    pass_co_exist = True
    pass_strand_bias = True
    pass_sequence_entropy = True
    match_count, ins_length = 0, 0
    hetero_germline_set = set()
    homo_germline_set = set()

    args.qual = args.qual if args.qual is not None else 102
    args.af = args.af if args.af is not None else 1.0

    disable_read_start_end_filtering = args.disable_read_start_end_filtering

    if not os.path.exists(tumor_bam_fn):
        tumor_bam_fn += ctg_name + '.bam'

    if args.hetero_info is not None and args.hetero_info != "":
        hetero_germline_set = set([tuple(item.split('-')) for item in args.hetero_info.split(',')])
    if args.homo_info is not None and args.homo_info != "":
        homo_germline_set = set([tuple(item.split('-')) for item in args.homo_info.split(',')])

    flanking = args.flanking

    ctg_range = "{}:{}-{}".format(ctg_name, max(pos - flanking, 1), pos + flanking + 1)
    samtools_command = "{} mpileup --min-MQ {} --min-BQ {} --excl-flags 2316 -r {} --output-MQ --output-QNAME --output-extra HP ".format(
        samtools, min_mq, min_bq, ctg_range)

    tumor_samtools_command = samtools_command + tumor_bam_fn

    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools,
        fasta_file_path=ref_fn,
        regions=[ctg_range]
    )

    # tumor
    pos_dict = defaultdict(defaultdict)
    pos_counter_dict = defaultdict(defaultdict)
    hap_dict = defaultdict(int)
    all_read_start_end_set = set()
    ALL_HAP_LIST = [0, 0, 0]
    HAP_LIST = [0, 0, 0]
    ALL_HAP_FORWARD_LIST = [0, 0, 0]
    ALL_HAP_REVERSE_LIST = [0, 0, 0]
    HAP_FORWARD_LIST = [0, 0, 0]
    HAP_REVERSE_LIST = [0, 0, 0]
    alt_base_read_name_set = set()
    homo_germline_pos_set = set([int(item[0]) for item in homo_germline_set])
    hetero_germline_pos_set = set([int(item[0]) for item in hetero_germline_set])

    samtools_mpileup_tumor_process = subprocess_popen(shlex.split(tumor_samtools_command), stderr=subprocess.PIPE)
    for row in samtools_mpileup_tumor_process.stdout:
        columns = row.split('\t')
        read_name_list = columns[7].split(',')

        base_counter, base_list, read_start_end_set = get_base_list(columns)

        for b_idx, base in enumerate(base_list):
            if base[0] == '#' or (base[0] >= 'a' and base[0] <= 'z'):
                read_name_list[b_idx] += '_1'  # reverse
            else:
                read_name_list[b_idx] += '_0'  # forward

        base_list = [[''.join(item[0]).upper()] + [item[1]] for item in base_list]

        p = int(columns[1])
        ctg = columns[0]
        ctg = p if args.ctg_name is not None else (ctg, p)
        if ctg in hetero_germline_pos_set or p == pos:
            phasing_info = columns[8].strip('\n').split(',')
            for hap_idx, hap in enumerate(phasing_info):
                if hap in '12' and read_name_list[hap_idx] not in hap_dict:
                    hap_dict[read_name_list[hap_idx]] = int(hap)

        # make union of all read start and end
        if len(read_start_end_set) >= len(base_list) * eps_rse:
            all_read_start_end_set = all_read_start_end_set.union(
                set([read_name_list[r_idx] for r_idx in read_start_end_set]))

        pos_dict[p] = dict(zip(read_name_list, base_list))
        center_ref_base = reference_sequence[p - max(pos - flanking, 1)]

        # discard low BQ & MQ variants
        average_min_bq = param.ont_min_bq
        average_min_mq = param.min_mq
        if p == pos:
            bq_list = [ord(qual) - 33 for qual in columns[5]]
            mq_list = [ord(qual) - 33 for qual in columns[6]]
            if is_snp:
                alt_base_bq_set = [bq for key, value, bq in zip(read_name_list, base_list, bq_list) if
                                   ''.join(value) == alt_base]
                alt_base_mq_set = [mq for key, value, mq in zip(read_name_list, base_list, mq_list) if
                                   ''.join(value) == alt_base]
            elif is_ins:
                alt_base_bq_set = [bq for key, value, bq in zip(read_name_list, base_list, bq_list) if
                                   ''.join(value).replace('+', '').upper() == alt_base and '+' in ''.join(value)]
                alt_base_mq_set = [mq for key, value, mq in zip(read_name_list, base_list, mq_list) if
                                   ''.join(value).replace('+', '').upper() == alt_base and '+' in ''.join(value)]
            elif is_del:
                alt_base_bq_set = [bq for key, value, bq in zip(read_name_list, base_list, bq_list) if
                                   len(ref_base) == len(value[1]) and '-' in value[1]]
                alt_base_mq_set = [mq for key, value, mq in zip(read_name_list, base_list, mq_list) if
                                   len(ref_base) == len(value[1]) and '-' in value[1]]

            if len(alt_base_bq_set) > 0 and sum(alt_base_bq_set) / len(alt_base_bq_set) <= average_min_bq:
                pass_bq = False

            if len(alt_base_mq_set) > 0 and sum(alt_base_mq_set) / len(alt_base_mq_set) <= average_min_mq:
                pass_mq = False

            for rn in read_name_list:
                ALL_HAP_LIST[hap_dict[rn]] += 1
                if rn.endswith('0'):
                    ALL_HAP_FORWARD_LIST[hap_dict[rn]] += 1
                elif rn.endswith('1'):
                    ALL_HAP_REVERSE_LIST[hap_dict[rn]] += 1

            if is_snp:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if ''.join(value) == alt_base])
            elif is_ins:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if ''.join(value).replace('+', '').upper() == alt_base and '+' in ''.join(value)])
            elif is_del:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if len(ref_base) == len(value[1]) and '-' in value[1]])

            for rn in alt_base_read_name_set:
                HAP_LIST[hap_dict[rn]] += 1
                if rn.endswith('0'):
                    HAP_FORWARD_LIST[hap_dict[rn]] += 1
                elif rn.endswith('1'):
                    HAP_REVERSE_LIST[hap_dict[rn]] += 1

        if len(base_counter) == 1 and base_counter[center_ref_base] > 0:
            continue
        pos_counter_dict[p] = base_counter

    samtools_mpileup_tumor_process.stdout.close()
    samtools_mpileup_tumor_process.wait()

    # near to read start end and have high overlap
    if not disable_read_start_end_filtering:
        if len(all_read_start_end_set.intersection(alt_base_read_name_set)) >= 0.3 * len(alt_base_read_name_set):
            pass_read_start_end = False

    alt_hap_counter = Counter([hap_dict[key] for key in alt_base_read_name_set])

    hp0, hp1, hp2 = alt_hap_counter[0], alt_hap_counter[1], alt_hap_counter[2]
    MAX = max(hp1, hp2)
    MIN = min(hp1, hp2)
    af = float(args.af)
    if is_snp and af < LOW_AF_SNV:
        if hp1 * hp2 > 0 and (MIN > args.min_alt_coverage or MAX / MIN <= 10):
            pass_hetero_both_side = False
    elif not is_snp and af < LOW_AF_INDEL:
        if hp1 * hp2 > 0 and (MIN > args.min_alt_coverage or MAX / MIN <= 10):
            pass_hetero_both_side = False

    is_phasable = hp1 * hp2 == 0 or (MAX / MIN >= 5 and (hp1 > args.min_alt_coverage or hp2 > args.min_alt_coverage))
    hap_index = 0 if not is_phasable else (1 if hp1 > hp2 else 2)

    # position with high overlap with current pos
    match_count = 0
    ins_length = 0
    all_base_dict = defaultdict(int)
    base_dict = defaultdict(int)
    alt_base_dict = defaultdict(int)

    for p, rb_dict in pos_dict.items():
        rb = reference_sequence[p - max(pos - flanking, 1)]
        read_alt_dict = pos_dict[p]
        if p == pos:
            continue

        ins_length += sum(
            [min(len(v[1]) - 1, flanking * 2) for key, v in read_alt_dict.items() if len(v[1]) > 3 and v[1][0] == '+'])

        for k, v in read_alt_dict.items():
            hap = hap_dict[k]
            base = ''.join(v)
            all_base_dict[hap] += 1
            if base.upper() != rb:
                base_dict[hap] += 1
                if k in alt_base_read_name_set:
                    alt_base_dict[hap] += 1

        read_alt_dict = pos_dict[p]
        inter_set = set(read_alt_dict.keys()).intersection(alt_base_read_name_set)
        alt_list = []

        for key in inter_set:
            base = ''.join(read_alt_dict[key]).upper()
            if base != rb and base not in '#*':
                alt_list.append(base)
        alt_base_counter = sorted(Counter(alt_list).items(), key=lambda x: x[1], reverse=True)

        upper_bound = 1 + eps
        lower_bound = 1 - eps
        if len(alt_list) == 0 or (alt_base_counter[0][1] >= len(alt_base_read_name_set) * upper_bound) or (
                alt_base_counter[0][1] <= len(alt_base_read_name_set) * lower_bound):
            continue

        # cal all alt count in current position
        if pos_counter_dict[p][alt_base_counter[0][0]] >= alt_base_counter[0][1] * upper_bound:
            continue

        match_count += 1

    if hap_index > 0:
        for p, ab in hetero_germline_set:
            p = int(p)
            rb = reference_sequence[p - pos + flanking]

            read_alt_dict = pos_dict[p]
            # snp
            if len(rb) == 1 and len(ab) == 1:
                overlap_count = set([key for key, value in read_alt_dict.items() if ''.join(value) == ab])
            # ins
            elif len(rb) == 1 and len(ab) > 1:
                overlap_count = set(
                    [key for key, value in read_alt_dict.items() if ab[:2] in value[1][1:] and len(value[1]) > 1])
            # del
            elif len(rb) > 1 and len(ab) == 1:
                overlap_count = set([key for key, value in read_alt_dict.items() if '-' in ''.join(value)])

            # in the same phased haplotype
            phased_overlap_set = set([key for key in overlap_count if hap_dict[key] == hap_index])
            if len(phased_overlap_set) == 0 or len(phased_overlap_set) * 2 < float(len(overlap_count)):
                continue

            inter_set = set([key for key in alt_base_read_name_set if hap_dict[key] == hap_index]).intersection(
                phased_overlap_set)
            # in the same phased haplotype while no overlapping
            if len(inter_set) == 0:
                pass_hetero = False
                break

    for p, ab in homo_germline_set:
        p = int(p)
        rb = reference_sequence[p - pos + flanking]
        read_alt_dict = pos_dict[p]

        # is the homo confident
        if len(rb) == 1 and len(ab) == 1:
            homo_alt_key = set([key for key, value in read_alt_dict.items() if ''.join(value) == ab])
        elif len(rb) == 1 and len(ab) > 1:
            homo_alt_key = set(
                [key for key, value in read_alt_dict.items() if (ab[1:2] in value[1][1:]) and len(value[1]) > 1])
        elif len(rb) > 1 and len(ab) == 1:
            homo_alt_key = set([key for key, value in read_alt_dict.items() if '-' in ''.join(value)])

        homo_hap_counter = Counter([hap_dict[key] for key in homo_alt_key])
        all_homo_hap_counter = Counter([hap_dict[key] for key in read_alt_dict.keys()])

        all_hp0, all_hp1, all_hp2 = all_homo_hap_counter[0], all_homo_hap_counter[1], all_homo_hap_counter[2]
        hp0, hp1, hp2 = homo_hap_counter[0], homo_hap_counter[1], homo_hap_counter[2]
        af = (hp0 + hp1 + hp2) / float(all_hp0 + all_hp1 + all_hp2) if (all_hp0 + all_hp1 + all_hp2) > 0 else 0.0

        def phasble(all_hap_list, hap_list):
            all_hp0, all_hp1, all_hp2 = all_hap_list
            if all_hp1 * all_hp2 == 0:
                return False
            hp0, hp1, hp2 = hap_list
            MAX = max(hp1, hp2)
            MIN = min(hp1, hp2)
            if hp1 * hp2 > 0 and MAX / MIN <= 10:
                return False
            return True

        is_homo_alt_phasable = phasble([all_hp0, all_hp1, all_hp2], [hp0, hp1, hp2])

        if af < min_hom_germline_af or is_homo_alt_phasable:
            continue

        inter_set = set(list(read_alt_dict.keys())).intersection(alt_base_read_name_set)
        if len(inter_set) == 0:
            continue
        if len(rb) == 1 and len(ab) == 1:
            overlap_count = set(
                [key for key, value in read_alt_dict.items() if ''.join(value) == ab and key in inter_set])
        elif len(rb) == 1 and len(ab) > 1:
            overlap_count = set([key for key, value in read_alt_dict.items() if
                                 (ab[1:2] in value[1][1:]) and len(value[1]) > 1 and key in inter_set])
        elif len(rb) > 1 and len(ab) == 1:
            overlap_count = set(
                [key for key, value in read_alt_dict.items() if '-' in ''.join(value) and key in inter_set])
        if len(overlap_count) == 0 or len(overlap_count) / len(inter_set) < eps:
            pass_homo = False
            break

    depth = sum(ALL_HAP_LIST) if sum(ALL_HAP_LIST) > 0 else 1

    if match_count >= max_co_exist_read_num or ins_length / depth > 3:
        pass_co_exist = False

    all_hp0, all_hp1, all_hp2 = ALL_HAP_LIST
    hp0, hp1, hp2 = HAP_LIST
    phaseable = all_hp1 * all_hp2 > 0 and hp1 * hp2 == 0 and (int(hp1) > args.min_alt_coverage or int(hp2) > args.min_alt_coverage)

    a0 = sum(HAP_FORWARD_LIST)
    a1 = sum(HAP_REVERSE_LIST)
    r0 = sum([x - y for x, y in zip(ALL_HAP_FORWARD_LIST, HAP_FORWARD_LIST)])
    r1 = sum([x - y for x, y in zip(ALL_HAP_REVERSE_LIST, HAP_REVERSE_LIST)])

    base_count_table = [[a0, r0], [a1, r1]]
    p_value = fisher_exact(base_count_table)
    if is_snp and p_value < 0.001 or (a0 == 0 or a1 == 0):
        pass_strand_bias = False
    elif not is_snp and p_value < 0.01 or (a0 == 0 or a1 == 0):
        pass_strand_bias = False

    if not is_snp:
        candidate_sequence_entropy = sqeuence_entropy_from(reference_sequence=reference_sequence)
        if candidate_sequence_entropy < sequence_entropy_threshold:
            pass_sequence_entropy = False

    pass_hap = (pass_hetero and pass_homo and pass_hetero_both_side and pass_read_start_end and pass_bq and pass_mq and pass_co_exist and pass_strand_bias and pass_sequence_entropy)

    print(' '.join([ctg_name, str(pos), str(pass_hap), str(phaseable), str(pass_hetero), str(pass_homo),
                    str(pass_read_start_end), str(pass_bq), str(pass_mq), str(pass_co_exist),
                    str(pass_hetero_both_side), str(pass_strand_bias), str(round(p_value, 5)),
                    str(pass_sequence_entropy)]))


def update_filter_info(args, key, row_str, phasable_set, fail_set, fail_pass_hetero_set, fail_pass_homo_set,
                       fail_pass_read_start_end_set, fail_pass_bq_set, fail_pass_mq_set, fail_pass_co_exist_set,
                       fail_pass_hetero_both_side_set, fail_pass_strand_bias_set, strand_bias_p_value_dict,
                       fail_pass_sequence_entropy_set):
    ctg_name = key[0] if args.ctg_name is None else args.ctg_name
    pos = key[1] if args.ctg_name is None else key
    k = (ctg_name, pos)
    columns = row_str.split('\t')

    is_candidate_filtered = 0
    phaseable = k in phasable_set

    if phaseable:
        columns[7] = 'H;' + columns[7]

    if k in fail_set:
        columns[5] = '0.0000'
        columns[6] = "LowQual"
        is_candidate_filtered = 1
    if k in fail_pass_bq_set:
        columns[6] += ";"
        columns[6] += "LowAltBQ"
    if k in fail_pass_mq_set:
        columns[6] += ";"
        columns[6] += "LowAltMQ"
    if k in fail_pass_read_start_end_set:
        columns[6] += ";"
        columns[6] += "ReadStartEnd"
    if k in fail_pass_co_exist_set:
        columns[6] += ";"
        columns[6] += "VariantCluster"
    if k in fail_pass_hetero_set or k in fail_pass_homo_set:
        columns[6] += ";"
        columns[6] += "NoAncestry"
    if k in fail_pass_hetero_both_side_set:
        columns[6] += ";"
        columns[6] += "MultiHap"
    if k in fail_pass_strand_bias_set:
        columns[6] += ";"
        columns[6] += "StrandBias"
    if k in fail_pass_sequence_entropy_set:
        columns[6] += ";"
        columns[6] += "LowSeqEntropy"
    if k in strand_bias_p_value_dict.keys():
        strand_bias_p_value = strand_bias_p_value_dict[k]
        columns[7] += ";"
        columns[7] += "SB={}".format(strand_bias_p_value)

    row_str = '\t'.join(columns)

    return row_str, is_candidate_filtered


def haplotype_filter(args):
    ctg_name = args.ctg_name
    threads = args.threads
    threads_low = max(1, int(threads * 4 / 5))
    apply_haplotype_filtering = args.apply_haplotype_filtering
    pileup_vcf_fn = args.pileup_vcf_fn
    germline_vcf_fn = args.germline_vcf_fn
    flanking = args.flanking
    output_dir = args.output_dir
    max_co_exist_read_num = args.min_alt_coverage
    is_indel = args.is_indel
    if not os.path.exists(output_dir):
        subprocess.run("mkdir -p {}".format(output_dir), shell=True)

    pileup_output_vcf_fn = args.output_vcf_fn
    if not apply_haplotype_filtering:
        subprocess.run("ln -sf {} {}".format(pileup_vcf_fn, pileup_output_vcf_fn), shell=True)
        return

    germine_input_vcf_reader = VcfReader(vcf_fn=germline_vcf_fn,
                                         ctg_name=ctg_name,
                                         show_ref=False,
                                         keep_row_str=False,
                                         filter_tag="PASS",
                                         save_header=False,
                                         skip_genotype=False)
    germine_input_vcf_reader.read_vcf()
    germline_input_variant_dict = germine_input_vcf_reader.variant_dict

    germline_gt_list = []
    for key in list(germline_input_variant_dict.keys()):
        if sum(germline_input_variant_dict[key].genotype) == 1:
            germline_gt_list.append((key, 1))
        elif sum(germline_input_variant_dict[key].genotype) == 2:
            germline_gt_list.append((key, 2))

    input_vcf_reader = VcfReader(vcf_fn=pileup_vcf_fn,
                                 ctg_name=ctg_name,
                                 show_ref=args.show_ref,
                                 keep_row_str=True,
                                 discard_indel=False if is_indel else True,
                                 filter_tag=args.input_filter_tag,
                                 save_header=True,
                                 keep_af=True)
    input_vcf_reader.read_vcf()
    pileup_variant_dict = input_vcf_reader.variant_dict

    input_variant_dict = defaultdict()
    for k, v in pileup_variant_dict.items():
        if v.filter != "PASS":
            continue
        if args.test_pos and k != args.test_pos:
            continue
        input_variant_dict[k] = v

    output_vcf_header = input_vcf_reader.header
    last_format_line = '##FORMAT=<ID=TU,Number=1,Type=Integer,Description="Count of T in the tumor BAM">'
    output_vcf_header = delete_lines_after(output_vcf_header, last_format_line)
    p_vcf_writer = VcfWriter(vcf_fn=pileup_output_vcf_fn,
                             ctg_name=ctg_name,
                             ref_fn=args.ref_fn,
                             header=output_vcf_header,
                             show_ref_calls=True)

    if not is_indel:
        hap_info_output_path = os.path.join(output_dir, "HAP_INFO_SNV")
    else:
        hap_info_output_path = os.path.join(output_dir, "HAP_INFO_INDEL")
    with open(hap_info_output_path, 'w') as f:
        for key, POS in input_variant_dict.items():
            ctg_name = args.ctg_name if args.ctg_name is not None else key[0]
            pos = key if args.ctg_name is not None else key[1]
            hetero_flanking_list = []
            homo_flanking_list = []
            for gk, gt in germline_gt_list:
                p = gk if args.ctg_name is not None else gk[1]
                if p > pos + flanking:
                    break
                if p > pos - flanking and p != pos:
                    alt_base = germline_input_variant_dict[gk].alternate_bases[0]
                    if gt == 1:
                        hetero_flanking_list.append('-'.join([str(p), str(alt_base)]))
                    else:
                        homo_flanking_list.append('-'.join([str(p), str(alt_base)]))
            POS.extra_infos = [set(hetero_flanking_list), set(homo_flanking_list)]
            info_list = [ctg_name, str(pos), POS.reference_bases, POS.alternate_bases[0], str(POS.af), str(POS.qual), \
                         ','.join(hetero_flanking_list), ','.join(homo_flanking_list)]
            f.write(' '.join(info_list) + '\n')

    file_directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    main_entry = os.path.join(file_directory, "clairs_to.py")

    parallel_command = "{} -C ' ' -j{} {} {} haplotype_filtering".format(args.parallel, threads_low, args.pypy3, main_entry)
    parallel_command += " --ctg_name {1}"
    parallel_command += " --pos {2}"
    parallel_command += " --ref_base {3}"
    parallel_command += " --alt_base {4}"
    parallel_command += " --af {5}"
    parallel_command += " --qual {6}"
    parallel_command += " --hetero_info {7}"
    parallel_command += " --homo_info {8}"
    parallel_command += " --samtools " + str(args.samtools)
    parallel_command += " --tumor_bam_fn " + str(args.tumor_bam_fn)
    parallel_command += " --ref_fn " + str(args.ref_fn)
    parallel_command += " --disable_read_start_end_filtering " + str(args.disable_read_start_end_filtering)
    parallel_command += " :::: " + str(hap_info_output_path)

    haplotype_filter_process = subprocess_popen(shlex.split(parallel_command))

    total_num = 0
    phasable_set = set()
    fail_set = set()
    fail_pass_hetero_set = set()
    fail_pass_homo_set = set()
    fail_pass_read_start_end_set = set()
    fail_pass_bq_set = set()
    fail_pass_mq_set = set()
    fail_pass_co_exist_set = set()
    fail_pass_hetero_both_side_set = set()
    fail_pass_strand_bias_set = set()
    fail_pass_sequence_entropy_set = set()
    strand_bias_p_value_dict = dict()

    for row in haplotype_filter_process.stdout:
        columns = row.rstrip().split()
        if len(columns) < 4:
            continue
        total_num += 1
        ctg_name, pos, pass_hap, phasable, pass_hetero, pass_homo, pass_read_start_end, pass_bq, pass_mq, pass_co_exist, pass_hetero_both_side, pass_strand_bias, strand_bias_p_value, pass_sequence_entropy = columns[:14]
        pos = int(pos)
        pass_hap = str2bool(pass_hap)
        pass_hetero = str2bool(pass_hetero)
        pass_homo = str2bool(pass_homo)
        pass_read_start_end = str2bool(pass_read_start_end)
        pass_bq = str2bool(pass_bq)
        pass_mq = str2bool(pass_mq)
        pass_co_exist = str2bool(pass_co_exist)
        pass_hetero_both_side = str2bool(pass_hetero_both_side)
        pass_strand_bias = str2bool(pass_strand_bias)
        pass_sequence_entropy = str2bool(pass_sequence_entropy)
        phasable = str2bool(phasable)
        if not pass_hap:
            fail_set.add((ctg_name, pos))
        if not pass_hetero:
            fail_pass_hetero_set.add((ctg_name, pos))
        if not pass_homo:
            fail_pass_homo_set.add((ctg_name, pos))
        if not pass_read_start_end:
            fail_pass_read_start_end_set.add((ctg_name, pos))
        if not pass_bq:
            fail_pass_bq_set.add((ctg_name, pos))
        if not pass_mq:
            fail_pass_mq_set.add((ctg_name, pos))
        if not pass_co_exist:
            fail_pass_co_exist_set.add((ctg_name, pos))
        if not pass_hetero_both_side:
            fail_pass_hetero_both_side_set.add((ctg_name, pos))
        if not pass_strand_bias:
            fail_pass_strand_bias_set.add((ctg_name, pos))
        if not pass_sequence_entropy:
            fail_pass_sequence_entropy_set.add((ctg_name, pos))
        strand_bias_p_value_dict[(ctg_name, pos)] = strand_bias_p_value
        if phasable:
            phasable_set.add((ctg_name, pos))
        if total_num > 0 and total_num % 1000 == 0:
            print("[INFO] Processing in {}, total processed positions: {}".format(ctg_name, total_num))

    haplotype_filter_process.stdout.close()
    haplotype_filter_process.wait()

    for key in sorted(pileup_variant_dict.keys()):
        row_str = pileup_variant_dict[key].row_str.rstrip()
        row_str, is_candidate_filtered = update_filter_info(args, key, row_str, phasable_set, fail_set,
                                                            fail_pass_hetero_set,
                                                            fail_pass_homo_set, fail_pass_read_start_end_set,
                                                            fail_pass_bq_set, fail_pass_mq_set, fail_pass_co_exist_set,
                                                            fail_pass_hetero_both_side_set, fail_pass_strand_bias_set,
                                                            strand_bias_p_value_dict, fail_pass_sequence_entropy_set)
        p_vcf_writer.vcf_writer.write(row_str + '\n')

    p_vcf_writer.close()

    print("[INFO] Total input calls: {}, filtered by all hard filers: {}".format(len(input_variant_dict), len(fail_set)))
    print("[INFO] Total input calls: {}, filtered by low alt bq: {}".format(len(input_variant_dict), len(fail_pass_bq_set)))
    print("[INFO] Total input calls: {}, filtered by low alt mq: {}".format(len(input_variant_dict), len(fail_pass_mq_set)))
    print("[INFO] Total input calls: {}, filtered by read start and end: {}".format(len(input_variant_dict), len(fail_pass_read_start_end_set)))
    print("[INFO] Total input calls: {}, filtered by variant cluster: {}".format(len(input_variant_dict), len(fail_pass_co_exist_set)))
    print("[INFO] Total input calls: {}, filtered by no ancestry: {}".format(len(input_variant_dict), len(fail_pass_hetero_set) + len(fail_pass_homo_set)))
    print("[INFO] Total input calls: {}, filtered by multi haplotypes: {}".format(len(input_variant_dict), len(fail_pass_hetero_both_side_set)))
    print("[INFO] Total input calls: {}, filtered by strand bias: {}".format(len(input_variant_dict), len(fail_pass_strand_bias_set)))
    print("[INFO] Total input calls: {}, filtered by low sequence entropy: {}".format(len(input_variant_dict), len(fail_pass_sequence_entropy_set)))


def main():
    parser = ArgumentParser(description="Haplotype filtering for long-read data")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="Output vcf file")

    parser.add_argument('--germline_vcf_fn', type=str, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output VCF directory")

    parser.add_argument('--hap_info_fn', type=str, default=None,
                        help="Hap Info File")

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
    # options for advanced users
    parser.add_argument('--apply_haplotype_filtering', type=str2bool, default=True,
                        help="EXPERIMENTAL: Apply haplotype filtering to the variant calls")

    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="Minimum number of reads supporting an alternative allele required for a somatic variant to be called. Default: %(default)d")

    parser.add_argument('--max_overlap_distance', type=int, default=100000,
                        help="EXPERIMENTAL: The largest window size for two somatic variants to be considered together for haplotype filtering. Default: %(default)d")

    ## filtering for Indel candidates
    parser.add_argument('--is_indel', action='store_true',
                        help=SUPPRESS)

    ## test using one position
    parser.add_argument('--test_pos', type=int, default=None,
                        help=SUPPRESS)

    ## flakning window size to process
    parser.add_argument('--flanking', type=int, default=100,
                        help=SUPPRESS)

    parser.add_argument('--add_phasing_info', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--is_happy_format', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--debug', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--hetero_info', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--homo_info', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--pos', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--ref_base', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--alt_base', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--af', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--qual', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--disable_read_start_end_filtering', type=str2bool, default=False,
                        help=SUPPRESS)

    global args
    args = parser.parse_args()

    if args.pos is None:
        haplotype_filter(args)
    else:
        haplotype_filter_per_pos(args)


if __name__ == "__main__":
    main()
