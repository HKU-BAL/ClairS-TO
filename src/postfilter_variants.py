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

HIGH_QUAL = 12
min_hom_germline_af = 0.75
eps = 0.5
eps_rse = 0.2

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


def postfilter_per_pos(args):
    pos = args.pos
    ctg_name = args.ctg_name
    ref_base = args.ref_base
    alt_base = args.alt_base
    tumor_bam_fn = args.tumor_bam_fn
    ref_fn = args.ref_fn
    samtools = args.samtools
    min_bq = args.min_bq
    min_mq = args.min_mq
    max_co_exist_read_num = args.min_alt_coverage
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    is_ins = len(ref_base) == 1 and len(alt_base) > 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    pass_read_start_end = True
    pass_co_exist = True
    pass_strand_bias = True

    args.qual = args.qual if args.qual is not None else 102
    args.af = args.af if args.af is not None else 1.0

    disable_read_start_end_filtering = args.disable_read_start_end_filtering

    if not os.path.exists(tumor_bam_fn):
        tumor_bam_fn += ctg_name + '.bam'

    flanking = args.flanking

    ctg_range = "{}:{}-{}".format(ctg_name, pos - flanking, pos + flanking + 1)
    samtools_command = "{} mpileup  --min-MQ {} --min-BQ {} --excl-flags 2316 -r {} --output-MQ --output-QNAME ".format(
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

        for idx in range(len(read_name_list)):
            hap_dict[read_name_list[idx]] = 0

        if len(read_start_end_set) >= len(base_list) * eps_rse:
            all_read_start_end_set = all_read_start_end_set.union(
                set([read_name_list[r_idx] for r_idx in read_start_end_set]))

        pos_dict[p] = dict(zip(read_name_list, base_list))
        center_ref_base = reference_sequence[p - pos + flanking]

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

    # near to read start end and have high overlap
    if not disable_read_start_end_filtering:
        if len(all_read_start_end_set.intersection(alt_base_read_name_set)) >= 0.3 * len(alt_base_read_name_set):
            pass_read_start_end = False

    match_count = 0
    ins_length = 0
    all_base_dict = defaultdict(int)
    base_dict = defaultdict(int)
    alt_base_dict = defaultdict(int)

    for p, rb_dict in pos_dict.items():
        rb = reference_sequence[p - pos + flanking]
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

    depth = sum(ALL_HAP_LIST) if sum(ALL_HAP_LIST) > 0 else 1
    if match_count >= max_co_exist_read_num or (ins_length / depth > 3 and float(args.qual) < HIGH_QUAL):
        pass_co_exist = False

    a0 = sum(HAP_FORWARD_LIST)
    a1 = sum(HAP_REVERSE_LIST)
    r0 = sum([x - y for x, y in zip(ALL_HAP_FORWARD_LIST, HAP_FORWARD_LIST)])
    r1 = sum([x - y for x, y in zip(ALL_HAP_REVERSE_LIST, HAP_REVERSE_LIST)])

    base_count_table = [[a0, r0], [a1, r1]]
    p_value = fisher_exact(base_count_table)
    if p_value < 0.001:
        pass_strand_bias = False

    pass_hard_filter = (pass_read_start_end and pass_co_exist and pass_strand_bias)
    
    print(' '.join([ctg_name, str(pos), str(pass_hard_filter), str(pass_read_start_end), str(pass_co_exist), str(pass_strand_bias), str(round(p_value, 5))]))


def update_filter_info(args, key, row_str, fail_pos_set, fail_pass_read_start_end_set, fail_pass_co_exist_set,
                       fail_pass_strand_bias_set, strand_bias_p_value_dict):
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
    if k in strand_bias_p_value_dict.keys():
        strand_bias_p_value = strand_bias_p_value_dict[k]
        columns[7] += ";"
        columns[7] += "SB={}".format(strand_bias_p_value)

    row_str = '\t'.join(columns)

    return row_str, is_candidate_filtered


def postfilter(args):
    ctg_name = args.ctg_name
    threads = args.threads
    threads_low = max(1, int(threads * 4 / 5))
    enable_postfilter = args.enable_postfilter
    pileup_vcf_fn = args.pileup_vcf_fn
    flanking = args.flanking
    output_dir = args.output_dir
    is_indel = args.is_indel
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
    with open(pf_info_output_path, 'w') as f:
        for key, POS in input_variant_dict.items():
            ctg_name = args.ctg_name if args.ctg_name is not None else key[0]
            pos = key if args.ctg_name is not None else key[1]
            info_list = [ctg_name, str(pos), POS.reference_bases, POS.alternate_bases[0], str(POS.af), str(POS.qual)]
            f.write(' '.join(info_list) + '\n')

    file_directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    main_entry = os.path.join(file_directory, "clairs_to.py")

    parallel_command = "{} -C ' ' -j{} {} {} postfilter_variants".format(args.parallel, threads_low, args.pypy3, main_entry)
    parallel_command += " --ctg_name {1}"
    parallel_command += " --pos {2}"
    parallel_command += " --ref_base {3}"
    parallel_command += " --alt_base {4}"
    parallel_command += " --af {5}"
    parallel_command += " --qual {6}"
    parallel_command += " --samtools " + str(args.samtools)
    parallel_command += " --tumor_bam_fn " + str(args.tumor_bam_fn)
    parallel_command += " --ref_fn " + str(args.ref_fn)
    parallel_command += " --disable_read_start_end_filtering " + str(args.disable_read_start_end_filtering)
    parallel_command += " :::: " + str(pf_info_output_path)

    postfilter_process = subprocess_popen(shlex.split(parallel_command))

    total_num = 0
    fail_set = set()
    fail_pass_read_start_end_set = set()
    fail_pass_co_exist_set = set()
    fail_pass_strand_bias_set = set()
    strand_bias_p_value_dict = dict()

    for row in postfilter_process.stdout:
        columns = row.rstrip().split()
        if len(columns) < 4:
            continue
        total_num += 1
        ctg_name, pos, pass_hard_filter, pass_read_start_end, pass_co_exist, pass_strand_bias, strand_bias_p_value = columns[:7]
        pos = int(pos)
        pass_hard_filter = str2bool(pass_hard_filter)
        pass_read_start_end = str2bool(pass_read_start_end)
        pass_co_exist = str2bool(pass_co_exist)
        pass_strand_bias = str2bool(pass_strand_bias)
        if not pass_hard_filter:
            fail_set.add((ctg_name, pos))
        if not pass_read_start_end:
            fail_pass_read_start_end_set.add((ctg_name, pos))
        if not pass_co_exist:
            fail_pass_co_exist_set.add((ctg_name, pos))
        if not pass_strand_bias:
            fail_pass_strand_bias_set.add((ctg_name, pos))
        strand_bias_p_value_dict[(ctg_name, pos)] = strand_bias_p_value
        if total_num > 0 and total_num % 1000 == 0:
            print("[INFO] Processing in {}, total processed positions: {}".format(ctg_name, total_num))

    postfilter_process.stdout.close()
    postfilter_process.wait()

    for key in sorted(pileup_variant_dict.keys()):
        row_str = pileup_variant_dict[key].row_str.rstrip()
        row_str, is_candidate_filtered = update_filter_info(args, key, row_str, fail_set, fail_pass_read_start_end_set,
                       fail_pass_co_exist_set, fail_pass_strand_bias_set, strand_bias_p_value_dict)
        p_vcf_writer.vcf_writer.write(row_str + '\n')

    p_vcf_writer.close()

    print("[INFO] Total input calls: {}, filtered by all hard filers: {}".format(len(input_variant_dict), len(fail_set)))
    print("[INFO] Total input calls: {}, filtered by read start and end: {}".format(len(input_variant_dict), len(fail_pass_read_start_end_set)))
    print("[INFO] Total input calls: {}, filtered by variant cluster: {}".format(len(input_variant_dict), len(fail_pass_co_exist_set)))
    print("[INFO] Total input calls: {}, filtered by strand bias: {}".format(len(input_variant_dict), len(fail_pass_strand_bias_set)))


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
    
    # options for advanced users
    parser.add_argument('--enable_postfilter', type=str2bool, default=True,
                        help="EXPERIMENTAL: Apply haplotype filtering to the variant calls")

    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="Minimum number of reads supporting an alternative allele required for a somatic variant to be called. Default: %(default)d")

    ## filtering for Indel candidates
    parser.add_argument('--is_indel', action='store_true',
                        help=SUPPRESS)

    ## test using one position
    parser.add_argument('--test_pos', type=int, default=None,
                        help=SUPPRESS)

    ## flakning window size to process
    parser.add_argument('--flanking', type=int, default=100,
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