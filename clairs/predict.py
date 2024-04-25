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
import numpy as np
import logging
import torch
import shlex

from time import time
from argparse import ArgumentParser, SUPPRESS
from threading import Thread
from sys import stderr
from subprocess import PIPE, run, Popen

sys.path.insert(0, '/autofs/bal13/lchen/ClairS-TO/ClairS-TO')

from clairs.call_variants import output_vcf_from_probability, OutputConfig
from clairs.model import CvT, BiGRU_NACGT, CvT_Indel, BiGRU_NACGT_Indel
from shared.utils import IUPAC_base_to_ACGT_base_dict as BASE2ACGT, BASIC_BASES, str2bool, file_path_from, log_error, \
    log_warning, subprocess_popen, TensorStdout
import shared.param as param


def batches_from(iterable, item_from, batch_size=1):
    iterable = iter(iterable)
    while True:
        chunk = []
        for _ in range(batch_size):
            try:
                chunk.append(item_from(next(iterable)))
            except StopIteration:
                yield chunk
                return
        yield chunk


def print_output_message(
        output_file,
        chromosome,
        position,
        reference_base,
        tumor_alt_info,
        input_forward_acgt_count_ori,
        input_reverse_acgt_count_ori,
        probabilities_a,
        probabilities_c,
        probabilities_g,
        probabilities_t,
        probabilities_i,
        probabilities_d,
        probabilities_na,
        probabilities_nc,
        probabilities_ng,
        probabilities_nt,
        probabilities_ni,
        probabilities_nd,
        extra_infomation_string="",
        disable_indel_calling=False
):
    global call_fn
    if call_fn is not None:
        output_vcf_from_probability(
            chromosome,
            position,
            reference_base,
            tumor_alt_info,
            input_forward_acgt_count_ori,
            input_reverse_acgt_count_ori,
            probabilities_a,
            probabilities_c,
            probabilities_g,
            probabilities_t,
            probabilities_i,
            probabilities_d,
            probabilities_na,
            probabilities_nc,
            probabilities_ng,
            probabilities_nt,
            probabilities_ni,
            probabilities_nd,
            output_config=output_config,
            vcf_writer=output_file,
            disable_indel_calling=disable_indel_calling
        )
    else:
        if not disable_indel_calling:
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                chromosome,
                position,
                reference_base,
                tumor_alt_info,
                input_forward_acgt_count_ori,
                input_reverse_acgt_count_ori,
                ' '.join(["{:0.8f}".format(x) for x in probabilities_a]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_c]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_g]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_t]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_i]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_d]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_na]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_nc]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_ng]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_nt]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_ni]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_nd]),
                extra_infomation_string
            ))
        else:
            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                chromosome,
                position,
                reference_base,
                tumor_alt_info,
                input_forward_acgt_count_ori,
                input_reverse_acgt_count_ori,
                ' '.join(["{:0.8f}".format(x) for x in probabilities_a]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_c]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_g]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_t]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_na]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_nc]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_ng]),
                ' '.join(["{:0.8f}".format(x) for x in probabilities_nt]),
                extra_infomation_string
            ))


def tensor_generator_from(tensor_file_path, batch_size, pileup=False, min_rescale_cov=None, platform='ont'):
    float_type = 'float32'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0
    if pileup:
        channel_size = param.pileup_channel_size
        tensor_shape = [param.no_of_positions, channel_size]
    else:
        tensor_shape = param.input_shape_dict[platform]
    prod_tensor_shape = np.prod(tensor_shape)

    def item_from(row):
        contig, coord, seq, tumor_tensor, tumor_alt_info, variant_type, ref_center= row.split("\t")[:7]
        ref_center = ref_center.strip()
        tumor_matrix = [float(item) for item in tumor_tensor.split()]

        if pileup:
            apply_normalize = False
            if min_rescale_cov is not None:
                tumor_coverage = float(tumor_alt_info.split('-')[0])
                tumor_rescale = float(min_rescale_cov) / tumor_coverage if tumor_coverage > min_rescale_cov else None
            else:
                tumor_rescale = None

            channel_size = param.pileup_channel_size
            tumor_channel_size = param.tumor_channel_size
            tensor = []
            for idx in range(param.no_of_positions):
                if apply_normalize:
                    tensor += [float(item) / tumor_coverage for item in
                               tumor_matrix[idx * channel_size: (idx + 1) * channel_size]]
                else:
                    if tumor_rescale is not None:
                        tensor += [item * tumor_rescale for item in
                                   tumor_matrix[idx * tumor_channel_size: (idx + 1) * tumor_channel_size]]
                    else:
                        tensor += tumor_matrix[idx * tumor_channel_size: (idx + 1) * tumor_channel_size]
        else:
            tumor_depth = len(tumor_matrix) // tensor_shape[1] // tensor_shape[2]
            tensor_depth = tumor_depth
            padding_depth = tensor_shape[0] - tensor_depth
            prefix_padding_depth = int(padding_depth / 2)
            suffix_padding_depth = padding_depth - int(padding_depth / 2)
            prefix_zero_padding = [0] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
            suffix_zero_padding = [0] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
            tensor = prefix_zero_padding + tumor_matrix + suffix_zero_padding
        tensor = np.array(tensor, dtype=np.dtype(float_type))

        pos = contig + ":" + coord + ":" + seq
        return tensor, pos, seq, tumor_alt_info, variant_type, ref_center

    for batch in batches_from(fo, item_from=item_from, batch_size=batch_size):
        tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
        positions = []
        tumor_alt_info_list = []
        variant_type_list = []
        ref_center_list = []
        for tensor, pos, seq, tumor_alt_info, variant_type, ref_center in batch:
            if seq[param.flankingBaseNum] not in "ACGT":
                continue
            tensors[len(positions)] = tensor
            positions.append(pos)
            tumor_alt_info_list.append(tumor_alt_info)
            variant_type_list.append(variant_type)
            ref_center_list.append(ref_center)

        current_batch_size = len(positions)
        X = np.reshape(tensors, ([batch_size] + tensor_shape))

        if processed_tensors > 0 and processed_tensors % 20000 == 0:
            print("Processed %d tensors" % processed_tensors, file=sys.stderr)

        processed_tensors += current_batch_size

        if current_batch_size <= 0:
            continue
        yield X[:current_batch_size], positions[:current_batch_size], tumor_alt_info_list[:current_batch_size], \
            variant_type_list[:current_batch_size], ref_center_list[:current_batch_size]

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()


def batch_output(output_file, batch_chr_pos_seq, tumor_alt_info_list, ref_center_list, input_list_forward_acgt_count_ori, input_list_reverse_acgt_count_ori, prediction_a, prediction_c,
                          prediction_g, prediction_t, prediction_i, prediction_d, prediction_na, prediction_nc, prediction_ng, prediction_nt, prediction_ni, prediction_nd, disable_indel_calling=False):
    batch_size = len(batch_chr_pos_seq)

    label_shape_cum = 2
    batch_probabilities_a = prediction_a[:, :label_shape_cum]
    batch_probabilities_c = prediction_c[:, :label_shape_cum]
    batch_probabilities_g = prediction_g[:, :label_shape_cum]
    batch_probabilities_t = prediction_t[:, :label_shape_cum]
    batch_probabilities_i = prediction_i[:, :label_shape_cum] if not disable_indel_calling else None
    batch_probabilities_d = prediction_d[:, :label_shape_cum] if not disable_indel_calling else None
    batch_probabilities_na = prediction_na[:, :label_shape_cum]
    batch_probabilities_nc = prediction_nc[:, :label_shape_cum]
    batch_probabilities_ng = prediction_ng[:, :label_shape_cum]
    batch_probabilities_nt = prediction_nt[:, :label_shape_cum]
    batch_probabilities_ni = prediction_ni[:, :label_shape_cum] if not disable_indel_calling else None
    batch_probabilities_nd = prediction_nd[:, :label_shape_cum] if not disable_indel_calling else None

    if len(batch_probabilities_a) != batch_size:
        sys.exit(
            "Inconsistent shape between input tensor and output predictions %d/%d" %
            (batch_size, len(batch_probabilities_a))
        )

    if disable_indel_calling:
        for (
                chr_pos_seq,
                tumor_alt_info,
                ref_center,
                input_forward_acgt_count_ori,
                input_reverse_acgt_count_ori,
                probabilities_a,
                probabilities_c,
                probabilities_g,
                probabilities_t,
                probabilities_na,
                probabilities_nc,
                probabilities_ng,
                probabilities_nt
        ) in zip(
            batch_chr_pos_seq,
            tumor_alt_info_list,
            ref_center_list,
            input_list_forward_acgt_count_ori,
            input_list_reverse_acgt_count_ori,
            batch_probabilities_a,
            batch_probabilities_c,
            batch_probabilities_g,
            batch_probabilities_t,
            batch_probabilities_na,
            batch_probabilities_nc,
            batch_probabilities_ng,
            batch_probabilities_nt
        ):
            probabilities_i = None
            probabilities_d = None
            probabilities_ni = None
            probabilities_nd = None
            output_with(
                output_file,
                chr_pos_seq,
                tumor_alt_info,
                ref_center,
                input_forward_acgt_count_ori,
                input_reverse_acgt_count_ori,
                probabilities_a,
                probabilities_c,
                probabilities_g,
                probabilities_t,
                probabilities_i,
                probabilities_d,
                probabilities_na,
                probabilities_nc,
                probabilities_ng,
                probabilities_nt,
                probabilities_ni,
                probabilities_nd,
                disable_indel_calling=disable_indel_calling
            )

    else:
        for (
                chr_pos_seq,
                tumor_alt_info,
                ref_center,
                input_forward_acgt_count_ori,
                input_reverse_acgt_count_ori,
                probabilities_a,
                probabilities_c,
                probabilities_g,
                probabilities_t,
                probabilities_i,
                probabilities_d,
                probabilities_na,
                probabilities_nc,
                probabilities_ng,
                probabilities_nt,
                probabilities_ni,
                probabilities_nd
        ) in zip(
            batch_chr_pos_seq,
            tumor_alt_info_list,
            ref_center_list,
            input_list_forward_acgt_count_ori,
            input_list_reverse_acgt_count_ori,
            batch_probabilities_a,
            batch_probabilities_c,
            batch_probabilities_g,
            batch_probabilities_t,
            batch_probabilities_i,
            batch_probabilities_d,
            batch_probabilities_na,
            batch_probabilities_nc,
            batch_probabilities_ng,
            batch_probabilities_nt,
            batch_probabilities_ni,
            batch_probabilities_nd
        ):
            output_with(
                output_file,
                chr_pos_seq,
                tumor_alt_info,
                ref_center,
                input_forward_acgt_count_ori,
                input_reverse_acgt_count_ori,
                probabilities_a,
                probabilities_c,
                probabilities_g,
                probabilities_t,
                probabilities_i,
                probabilities_d,
                probabilities_na,
                probabilities_nc,
                probabilities_ng,
                probabilities_nt,
                probabilities_ni,
                probabilities_nd,
                disable_indel_calling=disable_indel_calling
            )


def output_with(
        output_file,
        chr_pos_seq,
        tumor_alt_info,
        ref_center,
        input_forward_acgt_count_ori,
        input_reverse_acgt_count_ori,
        probabilities_a,
        probabilities_c,
        probabilities_g,
        probabilities_t,
        probabilities_i,
        probabilities_d,
        probabilities_na,
        probabilities_nc,
        probabilities_ng,
        probabilities_nt,
        probabilities_ni,
        probabilities_nd,
        disable_indel_calling=False
):
    if type(chr_pos_seq) == np.ndarray:
        chr_pos_seq = chr_pos_seq[0].decode()
        tumor_alt_info = tumor_alt_info[0].decode()
    elif type(chr_pos_seq) == np.bytes_ or type(chr_pos_seq) == bytes:
        chr_pos_seq = chr_pos_seq.decode()
        tumor_alt_info = tumor_alt_info.decode()

    chromosome, position, reference_sequence = chr_pos_seq.rstrip().split(':')[:3]
    reference_base = reference_sequence[param.flankingBaseNum].upper()
    print_output_message(
        output_file,
        chromosome,
        position,
        reference_base,
        tumor_alt_info,
        input_forward_acgt_count_ori,
        input_reverse_acgt_count_ori,
        probabilities_a,
        probabilities_c,
        probabilities_g,
        probabilities_t,
        probabilities_i,
        probabilities_d,
        probabilities_na,
        probabilities_nc,
        probabilities_ng,
        probabilities_nt,
        probabilities_ni,
        probabilities_nd,
        "",
        disable_indel_calling=disable_indel_calling
    )


def DataGenerator(dataset, num_epoch, batch_size, chunk_start_pos, chunk_end_pos):
    for idx in range(num_epoch):
        start_pos = chunk_start_pos + idx * batch_size
        end_pos = min(chunk_start_pos + (idx + 1) * batch_size, chunk_end_pos)
        input_matrix = dataset.input_matrix[start_pos:end_pos]
        position = dataset.position[start_pos:end_pos]  # .flatten()
        tumor_alt_info_list = dataset.tumor_alt_info[start_pos:end_pos]  # .flatten()
        ref_center_list =  dataset.ref_center[start_pos:end_pos]
        yield input_matrix, position, tumor_alt_info_list, ref_center_list


def predict(args):
    global output_config
    global call_fn

    output_config = OutputConfig(
        is_show_reference=args.show_ref,
        quality_score_for_pass=args.qual,
        pileup=args.pileup,
        disable_indel_calling=args.disable_indel_calling
    )

    param.flankingBaseNum = param.flankingBaseNum if args.flanking is None else args.flanking
    param.no_of_positions = param.flankingBaseNum * 2 + 1
    param.min_rescale_cov = param.min_rescale_cov if args.min_rescale_cov is None else args.min_rescale_cov
    predict_fn = args.predict_fn
    use_gpu = args.use_gpu
    variant_call_start_time = time()
    call_fn = args.call_fn
    chkpnt_fn_acgt = args.chkpnt_fn_acgt
    tensor_fn_acgt = args.tensor_fn_acgt
    chkpnt_fn_nacgt = args.chkpnt_fn_nacgt
    tensor_fn_nacgt = args.tensor_fn_nacgt
    platform = args.platform
    torch.set_num_threads(1)
    torch.manual_seed(0)
    np.random.seed(0)
    if use_gpu and not torch.cuda.is_available():
        print("[WARNING] --use_gpu is enabled, but cuda is not found")
        use_gpu = False
    if use_gpu:
        device = 'cuda'
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        device = 'cpu'

    if call_fn is not None:
        from shared.vcf import VcfWriter
        call_dir = os.path.dirname(call_fn)
        if not os.path.exists(call_dir):
            output = run("mkdir -p {}".format(call_dir), shell=True)
        vcf_writer = VcfWriter(vcf_fn=args.call_fn,
                               ref_fn=args.ref_fn,
                               show_ref_calls=args.show_ref,
                               sample_name=args.sample_name,
                               )
        output_file = vcf_writer
    elif predict_fn != "PIPE":
        predict_dir = os.path.dirname(predict_fn)
        if not os.path.exists(predict_dir):
            output = run("mkdir -p {}".format(predict_dir), shell=True)
        predict_fn_fpo = open(predict_fn, "wb")
        predict_fn_fp = subprocess_popen(shlex.split("{} -c".format(param.zstd)), stdin=PIPE, stdout=predict_fn_fpo)
        output_file = predict_fn_fp.stdin
    else:
        predict_fn_fp = TensorStdout(sys.stdout)
        output_file = predict_fn_fp.stdin

    global test_pos
    test_pos = None

    if args.disable_indel_calling:
        model_acgt = torch.load(chkpnt_fn_acgt, map_location=torch.device(device))
        model_aff = model_acgt['model_acgt']

        model_nacgt = torch.load(chkpnt_fn_nacgt, map_location=torch.device(device))
        model_neg = model_nacgt['model_nacgt']

    else:
        model_aff = CvT_Indel(
                num_classes=2,
                s1_emb_dim=16,  # stage 1 - dimension
                s1_emb_kernel=3,  # stage 1 - conv kernel
                s1_emb_stride=2,  # stage 1 - conv stride
                s1_proj_kernel=3,  # stage 1 - attention ds-conv kernel size
                s1_kv_proj_stride=2,  # stage 1 - attention key / value projection stride
                s1_heads=1,  # stage 1 - heads
                s1_depth=1,  # stage 1 - depth
                s1_mlp_mult=4,  # stage 1 - feedforward expansion factor
                s2_emb_dim=64,  # stage 2 - (same as above)
                s2_emb_kernel=3,
                s2_emb_stride=2,
                s2_proj_kernel=3,
                s2_kv_proj_stride=2,
                s2_heads=3,
                s2_depth=2,
                s2_mlp_mult=4,
                s3_emb_dim=128,  # stage 3 - (same as above)
                s3_emb_kernel=3,
                s3_emb_stride=2,
                s3_proj_kernel=3,
                s3_kv_proj_stride=2,
                s3_heads=4,
                s3_depth=3,
                s3_mlp_mult=4,
                dropout=0.,
                dropout_fc=0.3,
                depth=1,
                width=param.no_of_positions,
                dim=param.pileup_channel_size,
                apply_softmax=False,
                model_type="acgt"
            )

        model_acgt_saved = torch.load(chkpnt_fn_acgt, map_location=torch.device(device))
        model_aff_saved = model_acgt_saved['model_acgt']
        model_aff_state_dict = model_aff_saved.state_dict()
        model_aff.load_state_dict(model_aff_state_dict)

        model_neg = BiGRU_NACGT_Indel(apply_softmax=False,
                                 num_classes=2,
                                 channel_size=param.pileup_channel_size,
                                 model_type="nacgt")

        model_nacgt_saved = torch.load(chkpnt_fn_nacgt, map_location=torch.device(device))
        model_neg_saved = model_nacgt_saved['model_nacgt']
        model_neg_state_dict = model_neg_saved.state_dict()
        model_neg.load_state_dict(model_neg_state_dict)

    model_aff.eval()
    model_neg.eval()

    total = 0
    softmax = torch.nn.Softmax(dim=1)
    if not args.is_from_tables:
        is_finish_loaded_all_mini_batches = False
        mini_batches_loaded_acgt = []
        mini_batches_to_output_acgt = []
        mini_batches_loaded_nacgt = []
        mini_batches_to_output_nacgt = []
        mini_batches_loaded_acgt_ori = []
        mini_batches_to_output_acgt_ori = []

        def load_mini_batch():
            try:
                mini_batches_loaded_acgt.append(next(tensor_generator_acgt))
                mini_batches_loaded_nacgt.append(next(tensor_generator_nacgt))
                mini_batches_loaded_acgt_ori.append(next(tensor_generator_acgt_ori))
            except StopIteration:
                return

        tensor_generator_acgt = tensor_generator_from(tensor_file_path=tensor_fn_acgt,
                                                 batch_size=param.predictBatchSize,
                                                 pileup=args.pileup,
                                                 min_rescale_cov=param.min_rescale_cov,
                                                 platform=platform)

        tensor_generator_nacgt = tensor_generator_from(tensor_file_path=tensor_fn_nacgt,
                                                 batch_size=param.predictBatchSize,
                                                 pileup=args.pileup,
                                                 min_rescale_cov=param.min_rescale_cov,
                                                 platform=platform)

        tensor_generator_acgt_ori = tensor_generator_from(tensor_file_path=tensor_fn_acgt,
                                                 batch_size=param.predictBatchSize,
                                                 pileup=args.pileup,
                                                 min_rescale_cov=None,
                                                 platform=platform)

        while True:
            thread_pool = []
            if len(mini_batches_to_output_acgt) > 0:
                mini_batch_acgt = mini_batches_to_output_acgt.pop(0)
                input_tensor_acgt, position, tumor_alt_info_list, variant_type_list, ref_center_list = mini_batch_acgt

                mini_batch_nacgt = mini_batches_to_output_nacgt.pop(0)
                input_tensor_nacgt, position_nacgt, tumor_alt_info_list_nacgt, variant_type_list_nacgt, ref_center_list_nacgt = mini_batch_nacgt

                mini_batch_acgt_ori = mini_batches_to_output_acgt_ori.pop(0)
                input_tensor_acgt_ori, _, _, _, _ = mini_batch_acgt_ori

                input_matrix_acgt = torch.from_numpy(input_tensor_acgt).to(device)
                input_matrix_nacgt = torch.from_numpy(input_tensor_nacgt).to(device)
                input_matrix_acgt_ori = input_tensor_acgt_ori

                input_matrix_forward_acgt_count = input_matrix_acgt_ori[:, 16, 0:4]
                input_matrix_reverse_acgt_count = input_matrix_acgt_ori[:, 16, 9:13]
                input_matrix_forward_acgt_count_ori = input_matrix_forward_acgt_count.copy()
                input_matrix_reverse_acgt_count_ori = input_matrix_reverse_acgt_count.copy()

                negative_indices_forward = np.where(input_matrix_forward_acgt_count < 0)
                row_indices_forward, col_indices_forward = negative_indices_forward
                row_sums_forward = np.sum(input_matrix_forward_acgt_count[row_indices_forward], axis=1)
                input_matrix_forward_acgt_count_ori[row_indices_forward, col_indices_forward] = row_sums_forward * -1
                input_list_forward_acgt_count_ori = (np.where(input_matrix_forward_acgt_count_ori == -0, 0, input_matrix_forward_acgt_count_ori)).tolist()

                negative_indices_reverse = np.where(input_matrix_reverse_acgt_count < 0)
                row_indices_reverse, col_indices_reverse = negative_indices_reverse
                row_sums_reverse = np.sum(input_matrix_reverse_acgt_count[row_indices_reverse], axis=1)
                input_matrix_reverse_acgt_count_ori[row_indices_reverse, col_indices_reverse] = row_sums_reverse * -1
                input_list_reverse_acgt_count_ori = (np.where(input_matrix_reverse_acgt_count_ori == -0, 0,
                                                               input_matrix_reverse_acgt_count_ori)).tolist()

                with torch.no_grad():
                    if args.disable_indel_calling:
                        prediction_a, prediction_c, prediction_g, prediction_t = model_aff(
                            input_matrix_acgt)
                        prediction_na, prediction_nc, prediction_ng, prediction_nt = model_neg(
                            input_matrix_nacgt)
                        prediction_i = None
                        prediction_d = None
                        prediction_ni = None
                        prediction_nd = None
                    else:
                        prediction_a, prediction_c, prediction_g, prediction_t, prediction_i, prediction_d = model_aff(
                            input_matrix_acgt)
                        prediction_na, prediction_nc, prediction_ng, prediction_nt, prediction_ni, prediction_nd = model_neg(
                            input_matrix_nacgt)
                prediction_a = softmax(prediction_a)
                prediction_a = prediction_a.cpu().numpy()
                prediction_c = softmax(prediction_c)
                prediction_c = prediction_c.cpu().numpy()
                prediction_g = softmax(prediction_g)
                prediction_g = prediction_g.cpu().numpy()
                prediction_t = softmax(prediction_t)
                prediction_t = prediction_t.cpu().numpy()
                if not args.disable_indel_calling:
                    prediction_i = softmax(prediction_i)
                    prediction_i = prediction_i.cpu().numpy()
                    prediction_d = softmax(prediction_d)
                    prediction_d = prediction_d.cpu().numpy()
                prediction_na = softmax(prediction_na)
                prediction_na = prediction_na.cpu().numpy()
                prediction_nc = softmax(prediction_nc)
                prediction_nc = prediction_nc.cpu().numpy()
                prediction_ng = softmax(prediction_ng)
                prediction_ng = prediction_ng.cpu().numpy()
                prediction_nt = softmax(prediction_nt)
                prediction_nt = prediction_nt.cpu().numpy()
                if not args.disable_indel_calling:
                    prediction_ni = softmax(prediction_ni)
                    prediction_ni = prediction_ni.cpu().numpy()
                    prediction_nd = softmax(prediction_nd)
                    prediction_nd = prediction_nd.cpu().numpy()

                total += len(input_tensor_acgt)
                thread_pool.append(Thread(
                    target=batch_output,
                    args=(output_file, position, tumor_alt_info_list, ref_center_list, input_list_forward_acgt_count_ori, input_list_reverse_acgt_count_ori, prediction_a, prediction_c,
                          prediction_g, prediction_t, prediction_i, prediction_d, prediction_na, prediction_nc, prediction_ng, prediction_nt, prediction_ni, prediction_nd, args.disable_indel_calling)
                ))

            if not is_finish_loaded_all_mini_batches:
                thread_pool.append(Thread(target=load_mini_batch))

            for t in thread_pool:
                t.start()
            for t in thread_pool:
                t.join()

            is_finish_loaded_all_mini_batches = (len(mini_batches_loaded_acgt) == 0 and len(mini_batches_loaded_nacgt) == 0)
            while len(mini_batches_loaded_acgt) > 0:
                mini_batch_acgt = mini_batches_loaded_acgt.pop(0)
                mini_batches_to_output_acgt.append(mini_batch_acgt)

            while len(mini_batches_loaded_nacgt) > 0:
                mini_batch_nacgt = mini_batches_loaded_nacgt.pop(0)
                mini_batches_to_output_nacgt.append(mini_batch_nacgt)

            while len(mini_batches_loaded_acgt_ori) > 0:
                mini_batch_acgt_ori = mini_batches_loaded_acgt_ori.pop(0)
                mini_batches_to_output_acgt_ori.append(mini_batch_acgt_ori)

            is_nothing_to_predict_and_output = (
                    len(thread_pool) <= 0 and len(mini_batches_to_output_acgt) <= 0 and len(mini_batches_to_output_nacgt) <= 0
            )
            if is_finish_loaded_all_mini_batches and is_nothing_to_predict_and_output:
                break

    run_time = "%.1fs" % (time() - variant_call_start_time)
    logging.info("[INFO] {} total processed positions: {}, time elapsed: {}".format(args.ctg_name, total, run_time))

    if call_fn is not None:
        output_file.close()
        if os.path.exists(call_fn):
            vcf_file = open(call_fn, 'r').readlines()
            if not len(vcf_file):
                os.remove(call_fn)
            for row in vcf_file:
                if row[0] != '#':
                    return
            logging.info("[INFO] No variant output for {}, remove empty VCF".format(call_fn.split('/')[-1]))
            os.remove(call_fn)
    elif predict_fn != "PIPE":
        predict_fn_fp.stdin.close()
        predict_fn_fp.wait()
        predict_fn_fpo.close()


def main():
    parser = ArgumentParser(description="Candidate variants probability prediction using tensors and trained models")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--tensor_fn_acgt', type=str, default="PIPE",
                        help="Tensor input path of the affirmative model, or stdin if not set")

    parser.add_argument('--tensor_fn_nacgt', type=str, default="PIPE",
                        help="Tensor input path of the negational model, or stdin if not set")

    parser.add_argument('--chkpnt_fn_acgt', type=str, default=None,
                        help="Model path of the affirmative model, required")

    parser.add_argument('--chkpnt_fn_nacgt', type=str, default=None,
                        help="Model path of the negational model, required")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    # options for advanced users
    parser.add_argument('--min_rescale_cov', type=int, default=param.min_rescale_cov,
                        help="EXPERIMENTAL: Minimum coverage after rescalling from excessively high coverage data")

    parser.add_argument('--disable_indel_calling', type=str2bool, default=0,
                        help="EXPERIMENTAL: Disable Indel calling, default: enabled.")

    # options for debug purpose
    parser.add_argument('--predict_fn', type=str, default="PIPE",
                        help="DEBUG: Output network output probabilities for further analysis")

    # options for internal process control
    ## Use GPU for calling
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help=SUPPRESS)

    ## If set, variants with >=QUAL will be marked 'PASS', or 'LowQual'
    parser.add_argument('--qual', type=int, default=0,
                        help=SUPPRESS)

    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Use bin file from pytables to speed up calling.
    parser.add_argument('--is_from_tables', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--flanking', type=int, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    predict(args)


if __name__ == "__main__":
    main()
