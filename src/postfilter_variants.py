import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import shared.param as param
from shared.vcf import VcfReader, VcfWriter, Position
from shared.utils import str2bool, str_none, reference_sequence_from, subprocess_popen

import sys
import math
from shared.utils import IUPAC_base_to_num_dict as BASE2NUM

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
        if curP <= t:
            p_right_side += curP

    p += p_right_side

    return p


def calculate_sequence_entropy(sequence, entropy_window=None, kmer=5):

    count_repeat_kmer_counts = [0] * (entropy_window + 2)
    count_repeat_kmer_counts[0] = entropy_window

    entropy = [0.0] * (entropy_window + 2)
    for i in range(1, entropy_window + 2):
        e = 1.0 / entropy_window * i
        entropy[i] = e * math.log(e)
    entropy_mul = -1 / math.log(entropy_window)
    entropy_kmer_space = 1 << (2 * kmer)

    kmer_hash_counts = [0] * entropy_kmer_space
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
            kmer_prefix = ((kmer_prefix << 2) | n2) & mask
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
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':
            base_idx += 1
            read_start_set.add(len(base_list) - 1)
        if base == "$":
            read_end_set.add(len(base_list) - 1)
        base_idx += 1
    read_start_end_set = read_start_set if len(read_start_set) > len(read_end_set) else read_end_set
    upper_base_counter = Counter([''.join(item).upper() for item in base_list])
    return upper_base_counter, base_list, read_start_end_set


DEFAULT_POSTFILTER_CHUNK_MAX_SITES = 256
DEFAULT_POSTFILTER_CHUNK_MAX_SPAN = 50000
MIN_FANOUT_SITES_PER_MPILEUP_SHARD = 8


def _fanout_chunk_jobs_for_parallelism(chunk_jobs, threads_low, flanking):
    if threads_low <= 1 or not chunk_jobs:
        return chunk_jobs
    jobs = [(c, tuple(sorted(b)), lo, hi) for c, b, lo, hi in chunk_jobs]
    min_piece = MIN_FANOUT_SITES_PER_MPILEUP_SHARD
    safety = 0
    while len(jobs) < threads_low and safety < threads_low * 16:
        safety += 1
        idx = max(range(len(jobs)), key=lambda i: len(jobs[i][1]))
        contig, bs, _, _ = jobs[idx]
        n = len(bs)
        if n < 2 * min_piece:
            break
        mid = n // 2
        left, right = bs[:mid], bs[mid:]
        repl = []
        for sub in (left, right):
            slo = max(min(sub) - flanking, 1)
            shi = max(sub) + flanking + 1
            repl.append((contig, sub, slo, shi))
        jobs = jobs[:idx] + repl + jobs[idx + 1:]
    return jobs


def _group_postfilter_chunk_positions(sorted_positions, flanking, max_sites, max_span):
    groups = []
    n = len(sorted_positions)
    if n == 0:
        return groups
    max_sites = max(1, int(max_sites))
    max_span = max(1, int(max_span))
    i = 0
    while i < n:
        j = i
        lo = sorted_positions[i] - flanking
        hi = sorted_positions[i] + flanking
        while j + 1 < n:
            nlo = min(lo, sorted_positions[j + 1] - flanking)
            nhi = max(hi, sorted_positions[j + 1] + flanking)
            span = nhi - max(nlo, 1)
            if (j - i + 2) > max_sites or span > max_span:
                break
            lo, hi = nlo, nhi
            j += 1
        batch = sorted_positions[i : j + 1]
        reg_lo = max(lo, 1)
        reg_hi = hi + 1
        groups.append((batch, reg_lo, reg_hi))
        i = j + 1
    return groups


def _parse_mpileup_postfilter_chunk_dict(mpileup_stdout):
    chunk_rows = {}
    for row in mpileup_stdout:
        columns = row.split('\t')
        if len(columns) < 8:
            continue
        read_name_list = columns[7].split(',')
        base_counter, base_list, read_start_end_set = get_base_list(columns)
        for b_idx, base in enumerate(base_list):
            if base[0] == '#' or (base[0] >= 'a' and base[0] <= 'z'):
                read_name_list[b_idx] += '_1'
            else:
                read_name_list[b_idx] += '_0'
        base_list = [[''.join(item[0]).upper()] + [item[1]] for item in base_list]
        p = int(columns[1])
        chunk_rows[p] = {
            'read_name_list': read_name_list,
            'base_list': base_list,
            'base_counter': base_counter,
            'read_start_end_set': read_start_end_set,
        }
    return chunk_rows


def _run_mpileup_postfilter_chunk_dict(tumor_bam_fn, samtools, ctg_name, region_lo, region_hi, min_mq, min_bq):
    ctg_range = "{}:{}-{}".format(ctg_name, region_lo, region_hi)
    samtools_command = (
        "{} mpileup --min-MQ {} --min-BQ {} --excl-flags 2316 -r {} --output-MQ --output-QNAME ".format(
            samtools, min_mq, min_bq, ctg_range)
    )
    tumor_cmd = samtools_command + tumor_bam_fn
    proc = subprocess_popen(shlex.split(tumor_cmd), stderr=subprocess.PIPE)
    try:
        return _parse_mpileup_postfilter_chunk_dict(proc.stdout)
    finally:
        proc.stdout.close()
        proc.wait()


def _postfilter_finalize_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        pos_dict, pos_counter_dict, hap_dict,
        all_read_start_end_set, alt_base_read_name_set,
        ALL_HAP_LIST, HAP_LIST, ALL_HAP_FORWARD_LIST, ALL_HAP_REVERSE_LIST,
        HAP_FORWARD_LIST, HAP_REVERSE_LIST,
        ref_seq_site, ref_anchor,
        disable_read_start_end_filtering, max_co_exist_read_num):
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    pass_read_start_end = True
    pass_co_exist = True
    pass_strand_bias = True
    pass_sequence_entropy = True

    if not disable_read_start_end_filtering:
        if len(all_read_start_end_set.intersection(alt_base_read_name_set)) >= 0.3 * len(alt_base_read_name_set):
            pass_read_start_end = False

    match_count = 0
    ins_length = 0
    all_base_dict = defaultdict(int)
    base_dict = defaultdict(int)
    alt_base_dict = defaultdict(int)

    for p, rb_dict in pos_dict.items():
        rb = ref_seq_site[p - ref_anchor]
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

        if p not in pos_counter_dict:
            continue
        if pos_counter_dict[p][alt_base_counter[0][0]] >= alt_base_counter[0][1] * upper_bound:
            continue

        match_count += 1

    depth = sum(ALL_HAP_LIST) if sum(ALL_HAP_LIST) > 0 else 1
    if match_count >= max_co_exist_read_num or ins_length / depth > 3:
        pass_co_exist = False

    a0 = sum(HAP_FORWARD_LIST)
    a1 = sum(HAP_REVERSE_LIST)
    r0 = sum([x - y for x, y in zip(ALL_HAP_FORWARD_LIST, HAP_FORWARD_LIST)])
    r1 = sum([x - y for x, y in zip(ALL_HAP_REVERSE_LIST, HAP_REVERSE_LIST)])

    base_count_table = [[a0, r0], [a1, r1]]
    p_value = fisher_exact(base_count_table)
    if p_value < 0.001:
        pass_strand_bias = False

    if not is_snp:
        candidate_sequence_entropy = sqeuence_entropy_from(reference_sequence=ref_seq_site)
        if candidate_sequence_entropy < sequence_entropy_threshold:
            pass_sequence_entropy = False

    pass_hard_filter = (pass_read_start_end and pass_co_exist and pass_strand_bias and pass_sequence_entropy)

    return ' '.join([ctg_name, str(pos), str(pass_hard_filter), str(pass_read_start_end), str(pass_co_exist),
                     str(pass_strand_bias), str(round(p_value, 5)), str(pass_sequence_entropy)])


def _postfilter_build_state_and_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        chunk_rows, chunk_ref, region_lo,
        disable_read_start_end_filtering, max_co_exist_read_num):
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    is_ins = len(ref_base) == 1 and len(alt_base) > 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1

    ref_anchor = max(pos - flanking, 1)
    ref_end = pos + flanking + 1
    i0 = ref_anchor - region_lo
    i1 = ref_end - region_lo + 1
    ref_seq_site = chunk_ref[i0:i1]

    win_lo = max(pos - flanking, 1)
    win_hi = pos + flanking

    pos_dict = defaultdict(dict)
    pos_counter_dict = defaultdict(Counter)
    hap_dict = defaultdict(int)
    all_read_start_end_set = set()
    ALL_HAP_LIST = [0, 0, 0]
    HAP_LIST = [0, 0, 0]
    ALL_HAP_FORWARD_LIST = [0, 0, 0]
    ALL_HAP_REVERSE_LIST = [0, 0, 0]
    HAP_FORWARD_LIST = [0, 0, 0]
    HAP_REVERSE_LIST = [0, 0, 0]
    alt_base_read_name_set = set()

    for p in range(win_lo, win_hi + 1):
        if p not in chunk_rows:
            continue
        rec = chunk_rows[p]
        read_name_list = rec['read_name_list']
        base_list = rec['base_list']
        base_counter = rec['base_counter']
        read_start_end_set = rec['read_start_end_set']

        for idx in range(len(read_name_list)):
            hap_dict[read_name_list[idx]] = 0

        if len(read_start_end_set) >= len(base_list) * eps_rse:
            all_read_start_end_set = all_read_start_end_set.union(
                set([read_name_list[r_idx] for r_idx in read_start_end_set]))

        pos_dict[p] = dict(zip(read_name_list, base_list))
        center_ref_base = ref_seq_site[p - ref_anchor]

        if p == pos:
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
                    [key for key, value in zip(read_name_list, base_list) if
                     ''.join(value).replace('+', '').upper() == alt_base and '+' in ''.join(value)])
            elif is_del:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if
                     len(ref_base) == len(value[1]) and '-' in value[1]])

            for rn in alt_base_read_name_set:
                HAP_LIST[hap_dict[rn]] += 1
                if rn.endswith('0'):
                    HAP_FORWARD_LIST[hap_dict[rn]] += 1
                elif rn.endswith('1'):
                    HAP_REVERSE_LIST[hap_dict[rn]] += 1

        if len(base_counter) == 1 and base_counter[center_ref_base] > 0:
            continue
        pos_counter_dict[p] = base_counter

    pos_dict = dict(pos_dict)
    pos_counter_dict = dict(pos_counter_dict)

    return _postfilter_finalize_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        pos_dict, pos_counter_dict, hap_dict,
        all_read_start_end_set, alt_base_read_name_set,
        ALL_HAP_LIST, HAP_LIST, ALL_HAP_FORWARD_LIST, ALL_HAP_REVERSE_LIST,
        HAP_FORWARD_LIST, HAP_REVERSE_LIST,
        ref_seq_site, ref_anchor,
        disable_read_start_end_filtering, max_co_exist_read_num)


def postfilter_per_pos(args):
    pos = args.pos
    ctg_name = args.ctg_name
    ref_base = args.ref_base
    alt_base = args.alt_base
    tumor_bam_fn = args.tumor_bam_fn
    if not os.path.exists(tumor_bam_fn):
        tumor_bam_fn += ctg_name + '.bam'
    flanking = args.flanking
    region_lo = max(pos - flanking, 1)
    region_hi = pos + flanking + 1
    chunk_rows = _run_mpileup_postfilter_chunk_dict(
        tumor_bam_fn, args.samtools, ctg_name, region_lo, region_hi, args.min_mq, args.min_bq)
    chunk_ref = reference_sequence_from(
        samtools_execute_command=args.samtools,
        fasta_file_path=args.ref_fn,
        regions=["%s:%s-%s" % (ctg_name, region_lo, region_hi)])
    if chunk_ref is None:
        chunk_ref = ""
    line = _postfilter_build_state_and_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        chunk_rows, chunk_ref, region_lo,
        args.disable_read_start_end_filtering, args.min_alt_coverage)
    print(line)


def update_filter_info(args, key, row_str, fail_pos_set, fail_pass_read_start_end_set, fail_pass_co_exist_set,
                       fail_pass_strand_bias_set, strand_bias_p_value_dict, fail_pass_sequence_entropy_set):
    ctg_name = key[0] if args.ctg_name is None else args.ctg_name
    pos = key[1] if args.ctg_name is None else key
    k = (ctg_name, pos)
    columns = row_str.split('\t')

    is_candidate_filtered = 0

    if k in fail_pos_set:
        columns[5] = '0.0000'
        columns[6] = "LowQual"
        is_candidate_filtered = 1
    if k in fail_pass_read_start_end_set:
        columns[6] += ";"
        columns[6] += "ReadStartEnd"
    if k in fail_pass_co_exist_set:
        columns[6] += ";"
        columns[6] += "VariantCluster"
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


def _postfilter_candidates_progress(total_num):
    if total_num % 1000 == 0:
        print("[INFO] Postfilter variants: {} candidates processed".format(total_num), flush=True)


def _postfilter_per_pos_parallel_argv(
        args, threads_low, pf_info_path, main_entry, flanking, max_co_exist_read_num, is_indel):
    cmd_list = [
        args.parallel, '-C', ' ', '-j', str(threads_low),
        args.pypy3, main_entry, 'postfilter_variants',
        '--ctg_name', '{1}',
        '--pos', '{2}',
        '--ref_base', '{3}',
        '--alt_base', '{4}',
        '--af', '{5}',
        '--qual', '{6}',
        '--samtools', args.samtools,
        '--tumor_bam_fn', args.tumor_bam_fn,
        '--ref_fn', args.ref_fn,
        '--flanking', str(flanking),
        '--min_mq', str(args.min_mq),
        '--min_bq', str(args.min_bq),
        '--min_alt_coverage', str(max_co_exist_read_num),
        '--disable_read_start_end_filtering', str(args.disable_read_start_end_filtering),
        '--enable_postfilter', 'False',
        '::::',
        pf_info_path,
    ]
    if is_indel:
        idx = cmd_list.index('--enable_postfilter')
        cmd_list.insert(idx, '--is_indel')
    return cmd_list


def _ingest_postfilter_output_line(row, fail_set, fail_pass_read_start_end_set, fail_pass_co_exist_set,
                                   fail_pass_strand_bias_set, fail_pass_sequence_entropy_set,
                                   strand_bias_p_value_dict):
    columns = row.rstrip().split()
    if len(columns) < 8:
        return False
    ctg_o, pos, pass_hard_filter, pass_read_start_end, pass_co_exist, pass_strand_bias, strand_bias_p_value, pass_sequence_entropy = columns[:8]
    int_pos = int(pos)
    pass_hard_filter = str2bool(pass_hard_filter)
    pass_read_start_end = str2bool(pass_read_start_end)
    pass_co_exist = str2bool(pass_co_exist)
    pass_strand_bias = str2bool(pass_strand_bias)
    pass_sequence_entropy = str2bool(pass_sequence_entropy)
    if not pass_hard_filter:
        fail_set.add((ctg_o, int_pos))
    if not pass_read_start_end:
        fail_pass_read_start_end_set.add((ctg_o, int_pos))
    if not pass_co_exist:
        fail_pass_co_exist_set.add((ctg_o, int_pos))
    if not pass_strand_bias:
        fail_pass_strand_bias_set.add((ctg_o, int_pos))
    if not pass_sequence_entropy:
        fail_pass_sequence_entropy_set.add((ctg_o, int_pos))
    strand_bias_p_value_dict[(ctg_o, int_pos)] = strand_bias_p_value
    return True


def postfilter(args):
    ctg_name = args.ctg_name
    threads = args.threads
    threads_low = max(1, int(threads * 4 / 5))
    enable_postfilter = args.enable_postfilter
    pileup_vcf_fn = args.pileup_vcf_fn
    flanking = args.flanking
    output_dir = args.output_dir
    is_indel = args.is_indel
    max_co_exist_read_num = args.min_alt_coverage
    if not os.path.exists(output_dir):
        subprocess.run("mkdir -p {}".format(output_dir), shell=True)

    pileup_output_vcf_fn = args.output_vcf_fn
    if not enable_postfilter:
        subprocess.run("ln -sf {} {}".format(pileup_vcf_fn, pileup_output_vcf_fn), shell=True)
        return

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
        pf_info_output_path = os.path.join(output_dir, "PF_INFO_SNV")
    else:
        pf_info_output_path = os.path.join(output_dir, "PF_INFO_INDEL")
    _pf_buf_size = min(16 * 1024 * 1024, max(1024 * 1024, 1024 * (1 + len(input_variant_dict) // 10000)))
    sites_by_contig = defaultdict(list)
    site_meta = {}
    with open(pf_info_output_path, 'w', buffering=_pf_buf_size) as f:
        for key, POS in input_variant_dict.items():
            ctg_name_out = args.ctg_name if args.ctg_name is not None else key[0]
            pos = key if args.ctg_name is not None else key[1]
            site_meta[(ctg_name_out, pos)] = (POS.reference_bases, POS.alternate_bases[0])
            sites_by_contig[ctg_name_out].append(pos)
            info_list = [ctg_name_out, str(pos), POS.reference_bases, POS.alternate_bases[0], str(POS.af), str(POS.qual)]
            f.write(' '.join(info_list) + '\n')

    for _cname in sites_by_contig:
        sites_by_contig[_cname].sort()

    pass_site_count = sum(len(plist) for plist in sites_by_contig.values())
    variants_chunk_mode = getattr(args, 'postfilter_variants_chunk_mode', False)

    total_num = 0
    fail_set = set()
    fail_pass_read_start_end_set = set()
    fail_pass_co_exist_set = set()
    fail_pass_strand_bias_set = set()
    fail_pass_sequence_entropy_set = set()
    strand_bias_p_value_dict = dict()

    if variants_chunk_mode:
        chunk_max_sites = max(1, int(getattr(args, 'postfilter_chunk_max_sites', DEFAULT_POSTFILTER_CHUNK_MAX_SITES)))
        chunk_max_span = max(1, int(getattr(args, 'postfilter_chunk_max_span', DEFAULT_POSTFILTER_CHUNK_MAX_SPAN)))
        tumor_bam_work = args.tumor_bam_fn
        if ctg_name is not None and not os.path.exists(tumor_bam_work):
            tumor_bam_work = tumor_bam_work + ctg_name + '.bam'

        chunk_jobs = []
        for contig, plist in sites_by_contig.items():
            for batch, reg_lo, reg_hi in _group_postfilter_chunk_positions(
                    plist, flanking, chunk_max_sites, chunk_max_span):
                chunk_jobs.append((contig, tuple(batch), reg_lo, reg_hi))
        chunk_jobs = _fanout_chunk_jobs_for_parallelism(chunk_jobs, threads_low, flanking)

        def _run_postfilter_chunk(job):
            contig, batch, reg_lo, reg_hi = job
            chunk_rows = _run_mpileup_postfilter_chunk_dict(
                tumor_bam_work, args.samtools, contig, reg_lo, reg_hi, args.min_mq, args.min_bq)
            chunk_ref = reference_sequence_from(
                samtools_execute_command=args.samtools,
                fasta_file_path=args.ref_fn,
                regions=["{}:{}-{}".format(contig, reg_lo, reg_hi)])
            if chunk_ref is None:
                chunk_ref = ""
            out_lines = []
            for pos_p in batch:
                ref_b, alt_b = site_meta[(contig, pos_p)]
                out_lines.append(_postfilter_build_state_and_line(
                    contig, pos_p, ref_b, alt_b, flanking, chunk_rows, chunk_ref, reg_lo,
                    args.disable_read_start_end_filtering, max_co_exist_read_num))
            return out_lines

        if len(chunk_jobs) <= 1 or threads_low <= 1:
            for job in chunk_jobs:
                for line_str in _run_postfilter_chunk(job):
                    if _ingest_postfilter_output_line(
                            line_str, fail_set, fail_pass_read_start_end_set, fail_pass_co_exist_set,
                            fail_pass_strand_bias_set, fail_pass_sequence_entropy_set, strand_bias_p_value_dict):
                        total_num += 1
                        _postfilter_candidates_progress(total_num)
        else:
            with ThreadPoolExecutor(max_workers=min(threads_low, len(chunk_jobs))) as ex:
                for lines in ex.map(_run_postfilter_chunk, chunk_jobs):
                    for line_str in lines:
                        if _ingest_postfilter_output_line(
                                line_str, fail_set, fail_pass_read_start_end_set, fail_pass_co_exist_set,
                                fail_pass_strand_bias_set, fail_pass_sequence_entropy_set,
                                strand_bias_p_value_dict):
                            total_num += 1
                            _postfilter_candidates_progress(total_num)
    else:
        _src_dir = os.path.dirname(os.path.realpath(__file__))
        main_entry_pf = os.path.join(os.path.dirname(_src_dir), "clairs_to.py")
        if not os.path.isfile(main_entry_pf):
            print("[ERROR] clairs_to.py not found next to src/ at: {}".format(main_entry_pf), flush=True)
            sys.exit(1)
        if pass_site_count > 0:
            cmd_list = _postfilter_per_pos_parallel_argv(
                args, threads_low, pf_info_output_path, main_entry_pf, flanking,
                max_co_exist_read_num, is_indel)
            postfilter_parallel_proc = subprocess.Popen(
                cmd_list,
                stdout=subprocess.PIPE,
                stderr=None,
                universal_newlines=True,
                bufsize=8388608)
            for row in postfilter_parallel_proc.stdout:
                if _ingest_postfilter_output_line(
                        row, fail_set, fail_pass_read_start_end_set, fail_pass_co_exist_set,
                        fail_pass_strand_bias_set, fail_pass_sequence_entropy_set,
                        strand_bias_p_value_dict):
                    total_num += 1
                    _postfilter_candidates_progress(total_num)
            postfilter_parallel_proc.stdout.close()
            rc = postfilter_parallel_proc.wait()
            if rc != 0:
                print("[ERROR] postfilter variants per-position parallel (GNU parallel) exited with code {}".format(rc),
                      flush=True)
                sys.exit(rc)

    for key in sorted(pileup_variant_dict.keys()):
        row_str = pileup_variant_dict[key].row_str.rstrip()
        row_str, is_candidate_filtered = update_filter_info(args, key, row_str, fail_set,
                                                            fail_pass_read_start_end_set,
                                                            fail_pass_co_exist_set, fail_pass_strand_bias_set,
                                                            strand_bias_p_value_dict,
                                                            fail_pass_sequence_entropy_set)
        p_vcf_writer.vcf_writer.write(row_str + '\n')

    p_vcf_writer.close()


def main():
    parser = ArgumentParser(description="Post-filtering for short-read data")

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

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output VCF directory")

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

    parser.add_argument('--enable_postfilter', type=str2bool, default=True,
                        help="EXPERIMENTAL: Apply haplotype filtering to the variant calls")

    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="Minimum number of reads supporting an alternative allele required for a somatic variant to be called. Default: %(default)d")

    parser.add_argument('--is_indel', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--test_pos', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--flanking', type=int, default=100,
                        help=SUPPRESS)

    parser.add_argument('--postfilter_variants_chunk_mode', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--postfilter_chunk_max_sites', type=int, default=DEFAULT_POSTFILTER_CHUNK_MAX_SITES,
                        help=SUPPRESS)
    parser.add_argument('--postfilter_chunk_max_span', type=int, default=DEFAULT_POSTFILTER_CHUNK_MAX_SPAN,
                        help=SUPPRESS)

    parser.add_argument('--debug', action='store_true',
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
        postfilter(args)
    else:
        postfilter_per_pos(args)


if __name__ == "__main__":
    main()