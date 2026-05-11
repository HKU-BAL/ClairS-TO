import os
import shlex
import gc
import bisect
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import shared.param as param
from shared.vcf import VcfReader, VcfWriter, Position
from shared.utils import str2bool, str_none, reference_sequence_from

import sys
import math
from shared.utils import IUPAC_base_to_num_dict as BASE2NUM

LOW_AF_SNV = 0.1
LOW_AF_INDEL = 0.3
min_hom_germline_af = 0.75
# Match Clair-Mosaic ori_repo ClairS-TO/src/haplotype_filtering.py (variant cluster + homo overlap ratio).
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


DEFAULT_HAPLOTYPE_CHUNK_CANDIDATES_PER_MPILEUP = 200
DEFAULT_HAPLOTYPE_CHUNK_MAX_SPAN_BP = 5000000
MIN_FANOUT_SITES_PER_MPILEUP_SHARD = 8


def _partition_chunk_jobs_dual_cap(sites_by_contig, max_sites, max_span_bp, flanking):
    jobs = []
    ms = max(1, int(max_sites))
    mx_sp = int(max_span_bp)
    enforce_span = mx_sp > 0
    for contig, plist in sites_by_contig.items():
        sp = sorted(plist)
        n = len(sp)
        i = 0
        while i < n:
            j = i
            while j + 1 < n:
                trial = sp[i : j + 2]
                reg_lo = max(min(trial) - flanking, 1)
                reg_hi = max(trial) + flanking + 1
                span_bp = reg_hi - reg_lo
                if len(trial) > ms or (enforce_span and span_bp > mx_sp):
                    break
                j += 1
            batch = sp[i : j + 1]
            reg_lo = max(min(batch) - flanking, 1)
            reg_hi = max(batch) + flanking + 1
            jobs.append((contig, tuple(batch), reg_lo, reg_hi))
            i = j + 1
    return jobs


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


def _parse_mpileup_to_chunk_dict(mpileup_stdout):
    chunk_rows = {}
    for row in mpileup_stdout:
        columns = row.split('\t')
        if len(columns) < 9:
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
        phasing_parts = columns[8].strip('\n').split(',')
        bq_list = [ord(qual) - 33 for qual in columns[5]]
        mq_list = [ord(qual) - 33 for qual in columns[6]]
        chunk_rows[p] = {
            'read_name_list': read_name_list,
            'base_list': base_list,
            'base_counter': base_counter,
            'read_start_end_set': read_start_end_set,
            'phasing_parts': phasing_parts,
            'bq_list': bq_list,
            'mq_list': mq_list,
        }
    return chunk_rows


def _resolve_tumor_bam_path(tumor_bam_fn, chrom):
    if os.path.isfile(tumor_bam_fn):
        return tumor_bam_fn
    if tumor_bam_fn.endswith('.bam'):
        return tumor_bam_fn
    return tumor_bam_fn + chrom + '.bam'


def _merge_candidate_flank_bed_intervals(sorted_positions, flank):
    """UCSC BED: chromStart 0-based inclusive, chromEnd exclusive; covers 1-based genomic [win_lo, win_hi]."""
    merged = []
    for pos in sorted_positions:
        win_lo = max(pos - flank, 1)
        win_hi = pos + flank
        chrom_start = win_lo - 1
        chrom_end = win_hi
        if not merged:
            merged.append([chrom_start, chrom_end])
            continue
        ps, pe = merged[-1]
        if chrom_start <= pe:
            merged[-1][1] = max(pe, chrom_end)
        else:
            merged.append([chrom_start, chrom_end])
    return [(a[0], a[1]) for a in merged]


def _write_candidate_flank_bed(contig, positions, flank, bed_path):
    intervals = _merge_candidate_flank_bed_intervals(sorted(positions), flank)
    with open(bed_path, 'w') as bf:
        for s, e in intervals:
            bf.write("{}\t{}\t{}\n".format(contig, s, e))


def _run_mpileup_chunk_dict(tumor_bam_fn, samtools, ctg_name, region_lo, region_hi, min_mq, min_bq,
                           bed_file=None):
    tumor_cmd = (
        "{} mpileup --min-MQ {} --min-BQ {} --excl-flags 2316 ".format(
            shlex.quote(samtools), min_mq, min_bq))
    ctg_range = '{}:{}-{}'.format(ctg_name, region_lo, region_hi)
    if bed_file is not None:
        tumor_cmd += '-l {} '.format(shlex.quote(bed_file))
    tumor_cmd += '-r {} '.format(shlex.quote(ctg_range))
    tumor_cmd += '--output-MQ --output-QNAME --output-extra HP {}'.format(shlex.quote(tumor_bam_fn))

    proc = subprocess.Popen(
        shlex.split(tumor_cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        bufsize=8388608)
    try:
        chunk_rows = _parse_mpileup_to_chunk_dict(proc.stdout)
    finally:
        proc.stdout.close()
        mp_err = proc.stderr.read() or ''
        rc = proc.wait()
    if rc != 0:
        tail = mp_err.strip().replace('\n', ' ')
        if len(tail) > 800:
            tail = tail[:800] + '...'
        cmd_preview = tumor_cmd if len(tumor_cmd) <= 400 else tumor_cmd[:400] + '...'
        print(
            "[ERROR] samtools mpileup failed (exit {}). Command (trunc): {} stderr: {}".format(
                rc, cmd_preview, tail),
            flush=True)
    return chunk_rows


def _haplotype_finalize_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        pos_dict, pos_counter_dict, hap_dict,
        all_read_start_end_set, alt_base_read_name_set,
        ALL_HAP_LIST, HAP_LIST, ALL_HAP_FORWARD_LIST, ALL_HAP_REVERSE_LIST,
        HAP_FORWARD_LIST, HAP_REVERSE_LIST,
        hetero_germline_set, homo_germline_set,
        ref_seq_site, ref_anchor,
        disable_read_start_end_filtering, max_co_exist_read_num,
        af, qual,
        pass_bq, pass_mq):
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    is_ins = len(ref_base) == 1 and len(alt_base) > 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    pass_hetero = True
    pass_homo = True
    pass_hetero_both_side = True
    pass_read_start_end = True
    pass_co_exist = True
    pass_strand_bias = True
    pass_sequence_entropy = True
    match_count = 0
    ins_length = 0
    af = af if af is not None else 1.0

    if not disable_read_start_end_filtering:
        if len(alt_base_read_name_set) > 0:
            if len(all_read_start_end_set.intersection(alt_base_read_name_set)) >= 0.3 * len(
                    alt_base_read_name_set):
                pass_read_start_end = False

    alt_hap_counter = Counter([hap_dict[key] for key in alt_base_read_name_set])

    hp0, hp1, hp2 = alt_hap_counter[0], alt_hap_counter[1], alt_hap_counter[2]
    MAX = max(hp1, hp2)
    MIN = min(hp1, hp2)
    if is_snp and af < LOW_AF_SNV:
        if hp1 * hp2 > 0 and (MIN > max_co_exist_read_num or MAX / MIN <= 10):
            pass_hetero_both_side = False
    elif not is_snp and af < LOW_AF_INDEL:
        if hp1 * hp2 > 0 and (MIN > max_co_exist_read_num or MAX / MIN <= 10):
            pass_hetero_both_side = False

    is_phasable = hp1 * hp2 == 0 or (MAX / MIN >= 5 and (hp1 > max_co_exist_read_num or hp2 > max_co_exist_read_num))
    hap_index = 0 if not is_phasable else (1 if hp1 > hp2 else 2)

    all_base_dict = defaultdict(int)
    base_dict = defaultdict(int)
    alt_base_dict = defaultdict(int)

    for p, rb_dict in pos_dict.items():
        _ri = p - ref_anchor
        if _ri < 0 or _ri >= len(ref_seq_site):
            continue
        rb = ref_seq_site[_ri]
        read_alt_dict = pos_dict[p]
        # Match Clair-Mosaic ori_repo haplotype_filter_per_pos variant-cluster loop (only skip tumor pos).
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

    if hap_index > 0:
        for p, ab in hetero_germline_set:
            p = int(p)
            if p not in pos_dict:
                continue
            _ri = p - ref_anchor
            if _ri < 0 or _ri >= len(ref_seq_site):
                continue
            rb = ref_seq_site[_ri]

            read_alt_dict = pos_dict[p]
            if len(rb) == 1 and len(ab) == 1:
                overlap_count = set([key for key, value in read_alt_dict.items() if ''.join(value) == ab])
            elif len(rb) == 1 and len(ab) > 1:
                overlap_count = set(
                    [key for key, value in read_alt_dict.items() if ab[:2] in value[1][1:] and len(value[1]) > 1])
            elif len(rb) > 1 and len(ab) == 1:
                overlap_count = set([key for key, value in read_alt_dict.items() if '-' in ''.join(value)])
            else:
                overlap_count = set()

            phased_overlap_set = set([key for key in overlap_count if hap_dict[key] == hap_index])
            if len(phased_overlap_set) == 0 or len(phased_overlap_set) * 2 < float(len(overlap_count)):
                continue

            # Match Mosaic ori_repo: intersect suffixed pileup keys (not physical QNAME).
            inter_set = set(
                [key for key in alt_base_read_name_set if hap_dict[key] == hap_index]).intersection(
                phased_overlap_set)
            if len(inter_set) == 0:
                pass_hetero = False
                break

    for p, ab in homo_germline_set:
        p = int(p)
        if p not in pos_dict:
            continue
        _ri = p - ref_anchor
        if _ri < 0 or _ri >= len(ref_seq_site):
            continue
        rb = ref_seq_site[_ri]
        read_alt_dict = pos_dict[p]

        if len(rb) == 1 and len(ab) == 1:
            homo_alt_key = set([key for key, value in read_alt_dict.items() if ''.join(value) == ab])
        elif len(rb) == 1 and len(ab) > 1:
            homo_alt_key = set(
                [key for key, value in read_alt_dict.items() if (ab[1:2] in value[1][1:]) and len(value[1]) > 1])
        elif len(rb) > 1 and len(ab) == 1:
            homo_alt_key = set([key for key, value in read_alt_dict.items() if '-' in ''.join(value)])
        else:
            homo_alt_key = set()

        homo_hap_counter = Counter([hap_dict[key] for key in homo_alt_key])
        all_homo_hap_counter = Counter([hap_dict[key] for key in read_alt_dict.keys()])

        all_hp0, all_hp1, all_hp2 = all_homo_hap_counter[0], all_homo_hap_counter[1], all_homo_hap_counter[2]
        hp0, hp1, hp2 = homo_hap_counter[0], homo_hap_counter[1], homo_hap_counter[2]
        af_g = (hp0 + hp1 + hp2) / float(all_hp0 + all_hp1 + all_hp2) if (all_hp0 + all_hp1 + all_hp2) > 0 else 0.0

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

        if af_g < min_hom_germline_af or is_homo_alt_phasable:
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
        else:
            overlap_count = set()
        if len(overlap_count) == 0 or len(overlap_count) / len(inter_set) < eps:
            pass_homo = False
            break

    depth = sum(ALL_HAP_LIST) if sum(ALL_HAP_LIST) > 0 else 1

    if match_count >= max_co_exist_read_num or ins_length / depth > 3:
        pass_co_exist = False

    all_hp0, all_hp1, all_hp2 = ALL_HAP_LIST
    hp0, hp1, hp2 = HAP_LIST
    phaseable = all_hp1 * all_hp2 > 0 and hp1 * hp2 == 0 and (
            int(hp1) > max_co_exist_read_num or int(hp2) > max_co_exist_read_num)

    a0 = sum(HAP_FORWARD_LIST)
    a1 = sum(HAP_REVERSE_LIST)
    r0 = sum([x - y for x, y in zip(ALL_HAP_FORWARD_LIST, HAP_FORWARD_LIST)])
    r1 = sum([x - y for x, y in zip(ALL_HAP_REVERSE_LIST, HAP_REVERSE_LIST)])

    base_count_table = [[a0, r0], [a1, r1]]
    p_value = fisher_exact(base_count_table)
    # Match ori haplotype_filter_per_pos operator precedence (strict on empty alt strand counts).
    if is_snp and p_value < 0.001 or (a0 == 0 or a1 == 0):
        pass_strand_bias = False
    elif not is_snp and p_value < 0.01 or (a0 == 0 or a1 == 0):
        pass_strand_bias = False

    if not is_snp:
        candidate_sequence_entropy = sqeuence_entropy_from(reference_sequence=ref_seq_site)
        if candidate_sequence_entropy < sequence_entropy_threshold:
            pass_sequence_entropy = False

    pass_hap = (pass_hetero and pass_homo and pass_hetero_both_side and pass_read_start_end and pass_bq and pass_mq
                and pass_co_exist and pass_strand_bias and pass_sequence_entropy)

    return ' '.join([ctg_name, str(pos), str(pass_hap), str(phaseable), str(pass_hetero), str(pass_homo),
                     str(pass_read_start_end), str(pass_bq), str(pass_mq), str(pass_co_exist),
                     str(pass_hetero_both_side), str(pass_strand_bias), str(round(p_value, 5)),
                     str(pass_sequence_entropy)])


def _haplotype_build_state_and_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        chunk_rows, chunk_ref, region_lo,
        hetero_info_str, homo_info_str,
        disable_read_start_end_filtering, max_co_exist_read_num,
        af, qual):
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    is_ins = len(ref_base) == 1 and len(alt_base) > 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    ref_anchor = max(pos - flanking, 1)
    ref_end = pos + flanking + 1
    i0 = ref_anchor - region_lo
    i1 = ref_end - region_lo + 1
    if chunk_ref is None:
        chunk_ref = ''
    ref_seq_site = chunk_ref[i0:i1]

    hetero_germline_set = set()
    homo_germline_set = set()
    if hetero_info_str is not None and hetero_info_str != "":
        hetero_germline_set = set([tuple(item.split('-')) for item in hetero_info_str.split(',')])
    if homo_info_str is not None and homo_info_str != "":
        homo_germline_set = set([tuple(item.split('-')) for item in homo_info_str.split(',')])

    hetero_germline_pos_set = set([int(item[0]) for item in hetero_germline_set])

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
    pass_bq = True
    pass_mq = True

    for p in range(win_lo, win_hi + 1):
        if p not in chunk_rows:
            continue
        rec = chunk_rows[p]
        read_name_list = rec['read_name_list']
        base_list = rec['base_list']
        base_counter = rec['base_counter']
        read_start_end_set = rec['read_start_end_set']
        if p in hetero_germline_pos_set or p == pos:
            phasing_info = rec['phasing_parts']
            for hap_idx, hap in enumerate(phasing_info):
                if hap in '12' and read_name_list[hap_idx] not in hap_dict:
                    hap_dict[read_name_list[hap_idx]] = int(hap)

        if len(read_start_end_set) >= len(base_list) * eps_rse:
            all_read_start_end_set.update(read_name_list[r_idx] for r_idx in read_start_end_set)

        pos_dict[p] = dict(zip(read_name_list, base_list))

        if p == pos:
            bq_list = rec['bq_list']
            mq_list = rec['mq_list']
            average_min_bq = param.ont_min_bq
            average_min_mq = param.min_mq
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
            else:
                alt_base_bq_set = []
                alt_base_mq_set = []

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

        _ci = p - region_lo
        if not chunk_ref or _ci < 0 or _ci >= len(chunk_ref):
            continue
        center_ref_base = chunk_ref[_ci]
        if len(base_counter) == 1 and base_counter[center_ref_base] > 0:
            continue
        pos_counter_dict[p] = base_counter

    pos_dict = dict(pos_dict)
    pos_counter_dict = dict(pos_counter_dict)

    return _haplotype_finalize_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        pos_dict, pos_counter_dict, hap_dict,
        all_read_start_end_set, alt_base_read_name_set,
        ALL_HAP_LIST, HAP_LIST, ALL_HAP_FORWARD_LIST, ALL_HAP_REVERSE_LIST,
        HAP_FORWARD_LIST, HAP_REVERSE_LIST,
        hetero_germline_set, homo_germline_set,
        ref_seq_site, ref_anchor,
        disable_read_start_end_filtering, max_co_exist_read_num,
        af, qual, pass_bq, pass_mq)


def haplotype_filter_per_pos(args):
    pos = args.pos
    ctg_name = args.ctg_name
    if ctg_name is not None and ',' in ctg_name:
        print("[ERROR] haplotype_filter_per_pos requires a single chromosome in --ctg_name.", flush=True)
        sys.exit(1)
    ref_base = args.ref_base
    alt_base = args.alt_base
    tumor_bam_fn = args.tumor_bam_fn
    tumor_bam_fn = _resolve_tumor_bam_path(tumor_bam_fn, ctg_name)
    if not os.path.isfile(tumor_bam_fn):
        print("[ERROR] Tumor BAM not found: {}".format(tumor_bam_fn), flush=True)
        sys.exit(1)
    flanking = args.flanking
    region_lo = max(pos - flanking, 1)
    region_hi = pos + flanking + 1
    chunk_rows = _run_mpileup_chunk_dict(
        tumor_bam_fn, args.samtools, ctg_name, region_lo, region_hi, args.min_mq, args.min_bq,
        bed_file=None)
    chunk_ref = reference_sequence_from(
        samtools_execute_command=args.samtools,
        fasta_file_path=args.ref_fn,
        regions=["%s:%s-%s" % (ctg_name, region_lo, region_hi)])
    af = args.af if args.af is not None else 1.0
    qual = args.qual if args.qual is not None else 102
    line = _haplotype_build_state_and_line(
        ctg_name, pos, ref_base, alt_base, flanking,
        chunk_rows, chunk_ref, region_lo,
        args.hetero_info, args.homo_info,
        args.disable_read_start_end_filtering, args.min_alt_coverage, af, qual)
    print(line)


def update_filter_info(args, key, row_str, phasable_set, fail_set, fail_pass_hetero_set, fail_pass_homo_set,
                       fail_pass_read_start_end_set, fail_pass_bq_set, fail_pass_mq_set, fail_pass_co_exist_set,
                       fail_pass_hetero_both_side_set, fail_pass_strand_bias_set, strand_bias_p_value_dict,
                       fail_pass_sequence_entropy_set):
    if isinstance(key, tuple):
        ctg_name = key[0]
        pos = key[1]
    else:
        ctg_name = args.ctg_name
        pos = key
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


def _haplotype_candidates_progress(total_num):
    if total_num % 1000 == 0:
        print("[INFO] Haplotype filtering: {} candidates processed".format(total_num), flush=True)


def _haplotype_per_pos_parallel_argv(
        args, threads_low, hap_info_path, main_entry, flanking, max_co_exist_read_num, is_indel):
    cmd_list = [
        args.parallel, '-C', ' ', '-j', str(threads_low),
        args.pypy3, main_entry, 'haplotype_filtering',
        '--ctg_name', '{1}',
        '--pos', '{2}',
        '--ref_base', '{3}',
        '--alt_base', '{4}',
        '--af', '{5}',
        '--qual', '{6}',
        '--hetero_info', '{7}',
        '--homo_info', '{8}',
        '--samtools', args.samtools,
        '--tumor_bam_fn', args.tumor_bam_fn,
        '--ref_fn', args.ref_fn,
        '--flanking', str(flanking),
        '--min_mq', str(args.min_mq),
        '--min_bq', str(args.min_bq),
        '--min_alt_coverage', str(max_co_exist_read_num),
        '--disable_read_start_end_filtering', str(args.disable_read_start_end_filtering),
        '--apply_haplotype_filtering', 'False',
        '::::',
        hap_info_path,
    ]
    if is_indel:
        idx = cmd_list.index('--apply_haplotype_filtering')
        cmd_list.insert(idx, '--is_indel')
    return cmd_list


def _ingest_haplotype_filter_output_line(row, phasable_set, fail_set, fail_pass_hetero_set, fail_pass_homo_set,
                                         fail_pass_read_start_end_set, fail_pass_bq_set, fail_pass_mq_set,
                                         fail_pass_co_exist_set, fail_pass_hetero_both_side_set,
                                         fail_pass_strand_bias_set, fail_pass_sequence_entropy_set,
                                         strand_bias_p_value_dict):
    columns = row.rstrip().split()
    if len(columns) < 14:
        return False
    ctg_o, pos, pass_hap, phasable, pass_hetero, pass_homo, pass_read_start_end, pass_bq, pass_mq, pass_co_exist, pass_hetero_both_side, pass_strand_bias, strand_bias_p_value, pass_sequence_entropy = columns[:14]
    int_pos = int(pos)
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
        fail_set.add((ctg_o, int_pos))
    if not pass_hetero:
        fail_pass_hetero_set.add((ctg_o, int_pos))
    if not pass_homo:
        fail_pass_homo_set.add((ctg_o, int_pos))
    if not pass_read_start_end:
        fail_pass_read_start_end_set.add((ctg_o, int_pos))
    if not pass_bq:
        fail_pass_bq_set.add((ctg_o, int_pos))
    if not pass_mq:
        fail_pass_mq_set.add((ctg_o, int_pos))
    if not pass_co_exist:
        fail_pass_co_exist_set.add((ctg_o, int_pos))
    if not pass_hetero_both_side:
        fail_pass_hetero_both_side_set.add((ctg_o, int_pos))
    if not pass_strand_bias:
        fail_pass_strand_bias_set.add((ctg_o, int_pos))
    if not pass_sequence_entropy:
        fail_pass_sequence_entropy_set.add((ctg_o, int_pos))
    strand_bias_p_value_dict[(ctg_o, int_pos)] = strand_bias_p_value
    if phasable:
        phasable_set.add((ctg_o, int_pos))
    return True


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

    if ctg_name is not None and ',' not in ctg_name:
        _germline_entries = sorted(
            ((gk, gk, gt) for gk, gt in germline_gt_list),
            key=lambda t: t[0])
        _germline_pos_list = [t[0] for t in _germline_entries]
        _germline_by_contig = None
        _germline_pos_by_contig = None
    else:
        _buck = defaultdict(list)
        for gk, gt in germline_gt_list:
            if isinstance(gk, tuple):
                _chrom, _p = gk[0], gk[1]
            else:
                _chrom, _p = ctg_name, int(gk)
            _buck[_chrom].append((_p, gk, gt))
        _germline_by_contig = {}
        _germline_pos_by_contig = {}
        for _chrom, _lst in _buck.items():
            _lst.sort(key=lambda t: t[0])
            _germline_by_contig[_chrom] = _lst
            _germline_pos_by_contig[_chrom] = [t[0] for t in _lst]
        _germline_entries = None
        _germline_pos_list = None

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

    if args.input_filter_tag is not None:
        _hap_allowed_filters = frozenset(
            s.strip() for s in args.input_filter_tag.split(',') if s.strip())
    else:
        _hap_allowed_filters = None

    input_variant_dict = defaultdict()
    for k, v in pileup_variant_dict.items():
        if _hap_allowed_filters is not None:
            if v.filter not in _hap_allowed_filters:
                continue
        else:
            if v.filter != "PASS":
                continue
        if args.test_pos:
            _pk = k[1] if isinstance(k, tuple) else k
            if _pk != args.test_pos:
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
    _hap_buf_size = min(16 * 1024 * 1024, max(1024 * 1024, 1024 * (1 + len(input_variant_dict) // 10000)))
    sites_by_contig = defaultdict(list)
    site_meta = {}
    with open(hap_info_output_path, 'w', buffering=_hap_buf_size) as f:
        for key, POS in input_variant_dict.items():
            if isinstance(key, tuple):
                ctg_name_out, pos = key[0], key[1]
            else:
                ctg_name_out = args.ctg_name
                pos = key
            hetero_flanking_list = []
            homo_flanking_list = []
            if _germline_pos_list is not None:
                lo = bisect.bisect_right(_germline_pos_list, pos - flanking)
                hi = bisect.bisect_right(_germline_pos_list, pos + flanking)
                window = _germline_entries[lo:hi]
            else:
                _chrom = key[0]
                _pa = _germline_pos_by_contig.get(_chrom)
                if not _pa:
                    window = []
                else:
                    lo = bisect.bisect_right(_pa, pos - flanking)
                    hi = bisect.bisect_right(_pa, pos + flanking)
                    window = _germline_by_contig[_chrom][lo:hi]

            for p_gl, gk, gt in window:
                if p_gl == pos:
                    continue
                alt_base_g = germline_input_variant_dict[gk].alternate_bases[0]
                if gt == 1:
                    hetero_flanking_list.append('-'.join([str(p_gl), str(alt_base_g)]))
                else:
                    homo_flanking_list.append('-'.join([str(p_gl), str(alt_base_g)]))
            hetero_str = ','.join(hetero_flanking_list)
            homo_str = ','.join(homo_flanking_list)
            af_v = float(POS.af) if POS.af is not None else 1.0
            try:
                qual_v = float(POS.qual) if POS.qual is not None else 102.0
            except (TypeError, ValueError):
                qual_v = 102.0
            site_meta[(ctg_name_out, pos)] = (POS.reference_bases, POS.alternate_bases[0], af_v, qual_v, hetero_str, homo_str)
            sites_by_contig[ctg_name_out].append(pos)
            POS.extra_infos = [hetero_flanking_list, homo_flanking_list]
            info_list = [ctg_name_out, str(pos), POS.reference_bases, POS.alternate_bases[0], str(POS.af), str(POS.qual), \
                         hetero_str, homo_str]
            f.write(' '.join(info_list) + '\n')

    for _c in sites_by_contig:
        sites_by_contig[_c].sort()

    pass_site_count = sum(len(plist) for plist in sites_by_contig.values())
    chunk_mpileup = getattr(args, 'haplotype_filtering_chunk_mode', False)
    if sites_by_contig:
        for _c in sites_by_contig:
            _chrom = _c if (ctg_name is None or ',' in ctg_name) else ctg_name
            _bam_try = _resolve_tumor_bam_path(args.tumor_bam_fn, _chrom)
            if os.path.isfile(_bam_try):
                break
        else:
            _c0 = next(iter(sites_by_contig))
            _chrom0 = _c0 if (ctg_name is None or ',' in ctg_name) else ctg_name
            _p0 = _resolve_tumor_bam_path(args.tumor_bam_fn, _chrom0)
            print(
                "[ERROR] Tumor BAM not found. Phased prefix mode expects files like: prefix+<chrom>+.bam "
                "(run_clairs_to --phase_tumor). Example tried: {}".format(_p0),
                flush=True)
            sys.exit(1)

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

    if chunk_mpileup:
        chunk_max_sites = max(
            1,
            int(getattr(args, 'haplotype_chunk_max_sites', DEFAULT_HAPLOTYPE_CHUNK_CANDIDATES_PER_MPILEUP)))
        chunk_max_span_bp = int(getattr(args, 'haplotype_chunk_max_span', DEFAULT_HAPLOTYPE_CHUNK_MAX_SPAN_BP))

        chunk_jobs = _partition_chunk_jobs_dual_cap(
            sites_by_contig, chunk_max_sites, chunk_max_span_bp, flanking)
        chunk_jobs = _fanout_chunk_jobs_for_parallelism(chunk_jobs, threads_low, flanking)

        def _run_haplotype_chunk(job):
            contig, batch, _, _ = job
            bam_this = _resolve_tumor_bam_path(args.tumor_bam_fn, contig)
            batch_sorted = tuple(sorted(batch))
            sub_lo = max(min(batch_sorted) - flanking, 1)
            sub_hi = max(batch_sorted) + flanking + 1
            use_bed = getattr(args, 'haplotype_chunk_mpileup_bed', True)
            bed_path = None
            try:
                if use_bed:
                    safe_ctg = contig.replace('/', '_')
                    fd, bed_path = tempfile.mkstemp(suffix='.bed', prefix='hf_mpileup_{}_'.format(safe_ctg))
                    os.close(fd)
                    _write_candidate_flank_bed(contig, batch_sorted, flanking, bed_path)
                chunk_rows = _run_mpileup_chunk_dict(
                    bam_this, args.samtools, contig, sub_lo, sub_hi, args.min_mq, args.min_bq,
                    bed_file=bed_path)
                chunk_ref = reference_sequence_from(
                    samtools_execute_command=args.samtools,
                    fasta_file_path=args.ref_fn,
                    regions=["{}:{}-{}".format(contig, sub_lo, sub_hi)])
                if chunk_ref is None:
                    chunk_ref = ""
                out_lines = []
                for pos in batch_sorted:
                    ref_b, alt_b, af_v, qual_v, het_str, hom_str = site_meta[(contig, pos)]
                    out_lines.append(_haplotype_build_state_and_line(
                        contig, pos, ref_b, alt_b, flanking, chunk_rows, chunk_ref, sub_lo,
                        het_str, hom_str, args.disable_read_start_end_filtering, max_co_exist_read_num,
                        af_v, qual_v))
                return out_lines
            finally:
                if bed_path is not None:
                    try:
                        os.unlink(bed_path)
                    except OSError:
                        pass

        if len(chunk_jobs) <= 1 or threads_low <= 1:
            for job in chunk_jobs:
                lines = _run_haplotype_chunk(job)
                for line_str in lines:
                    if _ingest_haplotype_filter_output_line(
                            line_str, phasable_set, fail_set, fail_pass_hetero_set, fail_pass_homo_set,
                            fail_pass_read_start_end_set, fail_pass_bq_set, fail_pass_mq_set,
                            fail_pass_co_exist_set, fail_pass_hetero_both_side_set, fail_pass_strand_bias_set,
                            fail_pass_sequence_entropy_set, strand_bias_p_value_dict):
                        total_num += 1
                        _haplotype_candidates_progress(total_num)
        else:
            with ThreadPoolExecutor(max_workers=min(threads_low, len(chunk_jobs))) as ex:
                for lines in ex.map(_run_haplotype_chunk, chunk_jobs):
                    for line_str in lines:
                        if _ingest_haplotype_filter_output_line(
                                line_str, phasable_set, fail_set, fail_pass_hetero_set, fail_pass_homo_set,
                                fail_pass_read_start_end_set, fail_pass_bq_set, fail_pass_mq_set,
                                fail_pass_co_exist_set, fail_pass_hetero_both_side_set, fail_pass_strand_bias_set,
                                fail_pass_sequence_entropy_set, strand_bias_p_value_dict):
                            total_num += 1
                            _haplotype_candidates_progress(total_num)
    else:
        _src_dir = os.path.dirname(os.path.realpath(__file__))
        main_entry_hf = os.path.join(os.path.dirname(_src_dir), "clairs_to.py")
        if not os.path.isfile(main_entry_hf):
            print("[ERROR] clairs_to.py not found next to src/ at: {}".format(main_entry_hf), flush=True)
            sys.exit(1)
        if pass_site_count > 0:
            cmd_list = _haplotype_per_pos_parallel_argv(
                args, threads_low, hap_info_output_path, main_entry_hf, flanking,
                max_co_exist_read_num, is_indel)
            haplotype_filter_process = subprocess.Popen(
                cmd_list,
                stdout=subprocess.PIPE,
                stderr=None,
                universal_newlines=True,
                bufsize=8388608)
            for row in haplotype_filter_process.stdout:
                if _ingest_haplotype_filter_output_line(
                        row, phasable_set, fail_set, fail_pass_hetero_set, fail_pass_homo_set,
                        fail_pass_read_start_end_set, fail_pass_bq_set, fail_pass_mq_set,
                        fail_pass_co_exist_set, fail_pass_hetero_both_side_set, fail_pass_strand_bias_set,
                        fail_pass_sequence_entropy_set, strand_bias_p_value_dict):
                    total_num += 1
                    _haplotype_candidates_progress(total_num)
            haplotype_filter_process.stdout.close()
            rc = haplotype_filter_process.wait()
            if rc != 0:
                print("[ERROR] haplotype filtering per-position parallel (GNU parallel) exited with code {}".format(rc),
                      flush=True)
                sys.exit(rc)

    for key in sorted(pileup_variant_dict.keys()):
        row_str = pileup_variant_dict[key].row_str.rstrip()
        row_str, is_candidate_filtered = update_filter_info(args, key, row_str, phasable_set, fail_set,
                                                            fail_pass_hetero_set,
                                                            fail_pass_homo_set, fail_pass_read_start_end_set,
                                                            fail_pass_bq_set, fail_pass_mq_set, fail_pass_co_exist_set,
                                                            fail_pass_hetero_both_side_set, fail_pass_strand_bias_set,
                                                            strand_bias_p_value_dict, fail_pass_sequence_entropy_set)
        p_vcf_writer.vcf_writer.write(row_str + '\n')

    n_in = len(input_variant_dict)
    fail_no_anc = fail_pass_hetero_set | fail_pass_homo_set
    print("[INFO] Total input calls: {}, filtered by all hard filters: {}".format(n_in, len(fail_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by low alt bq: {}".format(n_in, len(fail_pass_bq_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by low alt mq: {}".format(n_in, len(fail_pass_mq_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by read start and end: {}".format(n_in, len(fail_pass_read_start_end_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by variant cluster: {}".format(n_in, len(fail_pass_co_exist_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by no ancestry: {}".format(
        n_in, len(fail_no_anc)), flush=True)
    print("[INFO] Total input calls: {}, filtered by multi haplotypes: {}".format(n_in, len(fail_pass_hetero_both_side_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by strand bias: {}".format(n_in, len(fail_pass_strand_bias_set)), flush=True)
    print("[INFO] Total input calls: {}, filtered by low sequence entropy: {}".format(n_in, len(fail_pass_sequence_entropy_set)), flush=True)

    p_vcf_writer.close()


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

    parser.add_argument('--is_indel', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--test_pos', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--flanking', type=int, default=100,
                        help=SUPPRESS)

    parser.add_argument('--haplotype_filtering_chunk_mode', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--haplotype_chunk_max_sites', type=int,
                        default=DEFAULT_HAPLOTYPE_CHUNK_CANDIDATES_PER_MPILEUP,
                        help=SUPPRESS)

    parser.add_argument('--haplotype_chunk_max_span', type=int,
                        default=DEFAULT_HAPLOTYPE_CHUNK_MAX_SPAN_BP,
                        help=SUPPRESS)

    parser.add_argument('--haplotype_chunk_mpileup_bed', type=str2bool, default=True,
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
