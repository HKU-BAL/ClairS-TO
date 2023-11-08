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

import torch
import numpy as np
from torch import nn, einsum
from einops import rearrange, repeat
import shared.param as param


def cal_scale(input_size, layers):
    output_size = input_size
    for _ in range(layers):
        output_size = np.ceil(output_size / 2)
    return int(output_size)

def group_dict_by_key(cond, d):
    return_val = [dict(), dict()]
    for key in d.keys():
        match = bool(cond(key))
        ind = int(not match)
        return_val[ind][key] = d[key]
    return (*return_val,)

def group_by_key_prefix_and_remove_prefix(prefix, d):
    kwargs_with_prefix, kwargs = group_dict_by_key(lambda x: x.startswith(prefix), d)
    kwargs_without_prefix = dict(map(lambda x: (x[0][len(prefix):], x[1]), tuple(kwargs_with_prefix.items())))
    return kwargs_without_prefix, kwargs

class LayerNorm(nn.Module): # layernorm, but done in the channel dimension #1
    def __init__(self, dim, eps = 1e-5):
        super().__init__()
        self.eps = eps
        self.g = nn.Parameter(torch.ones(1, dim, 1, 1))
        self.b = nn.Parameter(torch.zeros(1, dim, 1, 1))

    def forward(self, x):
        std = torch.var(x, dim = 1, unbiased = False, keepdim = True).sqrt()
        mean = torch.mean(x, dim = 1, keepdim = True)
        return (x - mean) / (std + self.eps) * self.g + self.b

class PreNorm(nn.Module):
    def __init__(self, dim, fn):
        super().__init__()
        self.norm = LayerNorm(dim)
        self.fn = fn
    def forward(self, x, **kwargs):
        x = self.norm(x)
        return self.fn(x, **kwargs)

class FeedForward(nn.Module):
    def __init__(self, dim, mult = 4, dropout = 0.):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(dim, dim * mult, 1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv2d(dim * mult, dim, 1),
            nn.Dropout(dropout)
        )
    def forward(self, x):
        return self.net(x)

class DepthWiseConv2d(nn.Module):
    def __init__(self, dim_in, dim_out, kernel_size, padding, stride, bias = True):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(dim_in, dim_in, kernel_size = kernel_size, padding = padding, groups = dim_in, stride = stride, bias = bias),
            nn.BatchNorm2d(dim_in),
            nn.Conv2d(dim_in, dim_out, kernel_size = 1, bias = bias)
        )
    def forward(self, x):
        return self.net(x)

class Attention(nn.Module):
    def __init__(self, dim, proj_kernel, kv_proj_stride, heads = 8, dim_head = 64, dropout = 0.):
        super().__init__()
        inner_dim = dim_head *  heads
        padding = proj_kernel // 2
        self.heads = heads
        self.scale = dim_head ** -0.5

        self.attend = nn.Softmax(dim = -1)

        self.to_q = DepthWiseConv2d(dim, inner_dim, proj_kernel, padding = padding, stride = 1, bias = False)
        self.to_kv = DepthWiseConv2d(dim, inner_dim * 2, proj_kernel, padding = padding, stride = kv_proj_stride, bias = False)

        self.to_out = nn.Sequential(
            nn.Conv2d(inner_dim, dim, 1),
            nn.Dropout(dropout)
        )

    def forward(self, x):
        shape = x.shape
        b, n, _, y, h = *shape, self.heads
        q, k, v = (self.to_q(x), *self.to_kv(x).chunk(2, dim = 1))
        q, k, v = map(lambda t: rearrange(t, 'b (h d) x y -> (b h) (x y) d', h = h), (q, k, v))

        dots = einsum('b i d, b j d -> b i j', q, k) * self.scale

        attn = self.attend(dots)

        out = einsum('b i j, b j d -> b i d', attn, v)
        out = rearrange(out, '(b h) (x y) d -> b (h d) x y', h = h, y = y)
        return self.to_out(out)

class Transformer(nn.Module):
    def __init__(self, dim, proj_kernel, kv_proj_stride, depth, heads, dim_head = 64, mlp_mult = 4, dropout = 0.):
        super().__init__()
        self.layers = nn.ModuleList([])
        for _ in range(depth):
            self.layers.append(nn.ModuleList([
                PreNorm(dim, Attention(dim, proj_kernel = proj_kernel, kv_proj_stride = kv_proj_stride, heads = heads, dim_head = dim_head, dropout = dropout)),
                PreNorm(dim, FeedForward(dim, mlp_mult, dropout = dropout))
            ]))
    def forward(self, x):
        for attn, ff in self.layers:
            x = attn(x) + x
            x = ff(x) + x
        return x

class CvT(nn.Module):
    def __init__(
        self,
        num_classes=2,
        s1_emb_dim=32,
        s1_emb_kernel=3,
        s1_emb_stride=2,
        s1_proj_kernel=3,
        s1_kv_proj_stride=2,
        s1_heads=1,
        s1_depth=1,
        s1_mlp_mult=4,
        s2_emb_dim=64,
        s2_emb_kernel=3,
        s2_emb_stride=2,
        s2_proj_kernel=3,
        s2_kv_proj_stride=2,
        s2_heads=3,
        s2_depth=2,
        s2_mlp_mult=4,
        s3_emb_dim=128,
        s3_emb_kernel=3,
        s3_emb_stride=2,
        s3_proj_kernel=3,
        s3_kv_proj_stride=2,
        s3_heads=6,
        s3_depth=10,
        s3_mlp_mult=4,
        dropout=0.,
        dropout_fc=0.3,
        depth=1,
        width=param.no_of_positions,
        dim=param.pileup_channel_size,
        apply_softmax=False,
        model_type="acgt"
    ):
        super().__init__()
        kwargs = dict(locals())
        dim = dim
        layers = []
        self.model_type = model_type
        self.layers_prefix = ('s1', 's2', 's3')
        for prefix in self.layers_prefix:
            config, kwargs = group_by_key_prefix_and_remove_prefix(f'{prefix}_', kwargs)
            layers.append(nn.Sequential(
                nn.Conv2d(dim, config['emb_dim'], kernel_size = config['emb_kernel'], padding = (config['emb_kernel'] // 2), stride = config['emb_stride']),
                LayerNorm(config['emb_dim']),
                Transformer(dim = config['emb_dim'], proj_kernel = config['proj_kernel'], kv_proj_stride = config['kv_proj_stride'], depth = config['depth'], heads = config['heads'], mlp_mult = config['mlp_mult'], dropout = dropout)
            ))

            dim = config['emb_dim']

        self.layer1 = layers[0]
        self.layer2 = layers[1]
        if len(self.layers_prefix) == 3:
            self.layer3 = layers[2]

        self.dropout_fc1 = nn.Dropout(p=dropout_fc)
        self.dropout_fc2 = nn.Dropout(p=dropout_fc)
        self.flatten = nn.Flatten()

        depth_scale_size = cal_scale(input_size=depth, layers=len(self.layers_prefix))
        width_scale_size = cal_scale(input_size=width, layers=len(self.layers_prefix))
        fc1_shape = dim * depth_scale_size * width_scale_size
        self.fc1 = nn.Linear(in_features=fc1_shape, out_features=128)
        self.fc2 = nn.Linear(in_features=128, out_features=num_classes)
        self.a_fc2 = nn.Linear(in_features=128, out_features=128)
        self.c_fc2 = nn.Linear(in_features=128, out_features=128)
        self.g_fc2 = nn.Linear(in_features=128, out_features=128)
        self.t_fc2 = nn.Linear(in_features=128, out_features=128)

        self.a_fc3 = nn.Linear(in_features=128, out_features=num_classes)
        self.c_fc3 = nn.Linear(in_features=128, out_features=num_classes)
        self.g_fc3 = nn.Linear(in_features=128, out_features=num_classes)
        self.t_fc3 = nn.Linear(in_features=128, out_features=num_classes)

        self.selu = nn.SELU()
        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)

    def forward(self, x):
        x = torch.unsqueeze(x, -2)
        x = x.permute(0, 3, 2, 1)
        output = self.layer1(x)
        output = self.layer2(output)
        if len(self.layers_prefix) == 3:
            output = self.layer3(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout_fc1(output)
        output = self.selu(self.dropout_fc1(self.fc1(output)))

        if self.model_type == "acgt":
            a_output = self.selu(self.dropout_fc2(self.a_fc2(output)))
            c_output = self.selu(self.dropout_fc2(self.c_fc2(output)))
            g_output = self.selu(self.dropout_fc2(self.g_fc2(output)))
            t_output = self.selu(self.dropout_fc2(self.t_fc2(output)))

            a_output = self.selu(self.a_fc3(a_output))
            c_output = self.selu(self.c_fc3(c_output))
            g_output = self.selu(self.g_fc3(g_output))
            t_output = self.selu(self.t_fc3(t_output))

            if self.apply_softmax:
                a_output = self.softmax(a_output)
                c_output = self.softmax(c_output)
                g_output = self.softmax(g_output)
                t_output = self.softmax(t_output)

            return a_output, c_output, g_output, t_output


class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, in_channels, out_channels, stride=1):
        super().__init__()

        # residual function
        self.residual_function = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels * BasicBlock.expansion, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(out_channels * BasicBlock.expansion)
        )

        # shortcut
        self.shortcut = nn.Sequential()

        if stride != 1 or in_channels != BasicBlock.expansion * out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels * BasicBlock.expansion, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels * BasicBlock.expansion)
            )

    def forward(self, x):
        return nn.ReLU(inplace=True)(self.residual_function(x) + self.shortcut(x))


class BiGRU(nn.Module):
    def __init__(
            self,
            num_classes=2,
            width=param.no_of_positions,
            batch_first=True,
            apply_softmax=False,
            channel_size=param.pileup_channel_size,
            model_type="acgt"
    ):
        super().__init__()
        kwargs = dict(locals())

        self.model_type = model_type

        self.num_layers = 2
        self.flatten = nn.Flatten()
        self.lstm_hidden_size = 128
        self.lstm_hidden_size2 = 192
        fc1_layer_size = 128
        dropout_rate = 0.3

        self.input_shape = [param.no_of_positions, self.lstm_hidden_size2 * 2]

        self.dim = channel_size
        self.lstm = nn.GRU(input_size=self.dim, hidden_size=self.lstm_hidden_size, batch_first=batch_first,
                           num_layers=1, bidirectional=True)

        self.lstm_2 = nn.GRU(input_size=self.lstm_hidden_size * 2, hidden_size=self.lstm_hidden_size2,
                             batch_first=batch_first,
                             num_layers=1, bidirectional=True)

        self.dropout_fc1 = nn.Dropout(p=dropout_rate)
        self.fc1 = nn.Linear(in_features=self.input_shape[0] * self.input_shape[1], out_features=fc1_layer_size)
        self.dropout_fc2 = nn.Dropout(p=dropout_rate)

        self.fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)

        self.a_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.c_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.g_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.t_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)

        self.na_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.nc_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.ng_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.nt_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)

        self.selu = nn.SELU()
        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)

    def forward(self, x):

        output, hidden = self.lstm(x)
        output, hidden = self.lstm_2(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout_fc1(output)
        output = self.selu(self.dropout_fc1(self.fc1(output)))

        output = self.selu(self.dropout_fc2(self.fc2(output)))

        if self.model_type == "acgt":

            a_output = self.selu(self.a_fc3(output))
            c_output = self.selu(self.c_fc3(output))
            g_output = self.selu(self.g_fc3(output))
            t_output = self.selu(self.t_fc3(output))

            return a_output, c_output, g_output, t_output

        elif self.model_type == "nacgt":

            na_output = self.selu(self.na_fc3(output))
            nc_output = self.selu(self.nc_fc3(output))
            ng_output = self.selu(self.ng_fc3(output))
            nt_output = self.selu(self.nt_fc3(output))

            if self.apply_softmax:
                na_output = self.softmax(na_output)
                nc_output = self.softmax(nc_output)
                ng_output = self.softmax(ng_output)
                nt_output = self.softmax(nt_output)

            return na_output, nc_output, ng_output, nt_output


class BiGRU_ACGT(nn.Module):
    def __init__(
            self,
            num_classes=2,
            width=param.no_of_positions,
            batch_first=True,
            apply_softmax=False,
            channel_size=param.pileup_channel_size,
            model_type="acgt"
    ):
        super().__init__()
        kwargs = dict(locals())

        self.model_type = model_type

        self.num_layers = 2
        self.flatten = nn.Flatten()
        self.lstm_hidden_size = 128
        self.lstm_hidden_size2 = 192
        fc1_layer_size = 128
        dropout_rate = 0.3

        self.input_shape = [param.no_of_positions, self.lstm_hidden_size2 * 2]

        self.dim = channel_size
        self.lstm = nn.GRU(input_size=self.dim, hidden_size=self.lstm_hidden_size, batch_first=batch_first,
                           num_layers=1, bidirectional=True)

        self.lstm_2 = nn.GRU(input_size=self.lstm_hidden_size * 2, hidden_size=self.lstm_hidden_size2,
                             batch_first=batch_first,
                             num_layers=1, bidirectional=True)

        self.dropout_fc1 = nn.Dropout(p=dropout_rate)
        self.fc1 = nn.Linear(in_features=self.input_shape[0] * self.input_shape[1], out_features=fc1_layer_size)
        self.dropout_fc2 = nn.Dropout(p=dropout_rate)

        self.fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)

        self.a_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)
        self.c_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)
        self.g_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)
        self.t_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)

        self.a_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.c_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.g_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.t_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)

        self.selu = nn.SELU()
        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)

    def forward(self, x):

        output, hidden = self.lstm(x)
        output, hidden = self.lstm_2(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout_fc1(output)
        output = self.selu(self.dropout_fc1(self.fc1(output)))

        if self.model_type == "acgt":
            a_output = self.selu(self.dropout_fc2(self.a_fc2(output)))
            c_output = self.selu(self.dropout_fc2(self.c_fc2(output)))
            g_output = self.selu(self.dropout_fc2(self.g_fc2(output)))
            t_output = self.selu(self.dropout_fc2(self.t_fc2(output)))

            a_output = self.selu(self.a_fc3(a_output))
            c_output = self.selu(self.c_fc3(c_output))
            g_output = self.selu(self.g_fc3(g_output))
            t_output = self.selu(self.t_fc3(t_output))

            if self.apply_softmax:
                a_output = self.softmax(a_output)
                c_output = self.softmax(c_output)
                g_output = self.softmax(g_output)
                t_output = self.softmax(t_output)

            return a_output, c_output, g_output, t_output


class BiGRU_NACGT(nn.Module):
    def __init__(
            self,
            num_classes=2,
            width=param.no_of_positions,
            batch_first=True,
            apply_softmax=False,
            channel_size=param.pileup_channel_size,
            model_type="acgt"
    ):
        super().__init__()
        kwargs = dict(locals())

        self.model_type = model_type

        self.num_layers = 2
        self.flatten = nn.Flatten()
        self.lstm_hidden_size = 128
        self.lstm_hidden_size2 = 192
        fc1_layer_size = 128
        dropout_rate = 0.3

        self.input_shape = [param.no_of_positions, self.lstm_hidden_size2 * 2]

        self.dim = channel_size
        self.lstm = nn.GRU(input_size=self.dim, hidden_size=self.lstm_hidden_size, batch_first=batch_first,
                           num_layers=1, bidirectional=True)

        self.lstm_2 = nn.GRU(input_size=self.lstm_hidden_size * 2, hidden_size=self.lstm_hidden_size2,
                             batch_first=batch_first,
                             num_layers=1, bidirectional=True)

        self.dropout_fc1 = nn.Dropout(p=dropout_rate)
        self.fc1 = nn.Linear(in_features=self.input_shape[0] * self.input_shape[1], out_features=fc1_layer_size)
        self.dropout_fc2 = nn.Dropout(p=dropout_rate)

        self.fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)

        self.na_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)
        self.nc_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)
        self.ng_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)
        self.nt_fc2 = nn.Linear(in_features=fc1_layer_size, out_features=fc1_layer_size)

        self.na_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.nc_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.ng_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)
        self.nt_fc3 = nn.Linear(in_features=fc1_layer_size, out_features=num_classes)

        self.selu = nn.SELU()
        self.apply_softmax = apply_softmax
        if self.apply_softmax:
            self.softmax = nn.Softmax(dim=-1)

    def forward(self, x):

        output, hidden = self.lstm(x)
        output, hidden = self.lstm_2(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout_fc1(output)
        output = self.selu(self.dropout_fc1(self.fc1(output)))

        if self.model_type == "nacgt":
            na_output = self.selu(self.dropout_fc2(self.na_fc2(output)))
            nc_output = self.selu(self.dropout_fc2(self.nc_fc2(output)))
            ng_output = self.selu(self.dropout_fc2(self.ng_fc2(output)))
            nt_output = self.selu(self.dropout_fc2(self.nt_fc2(output)))

            na_output = self.selu(self.na_fc3(na_output))
            nc_output = self.selu(self.nc_fc3(nc_output))
            ng_output = self.selu(self.ng_fc3(ng_output))
            nt_output = self.selu(self.nt_fc3(nt_output))

            if self.apply_softmax:
                na_output = self.softmax(na_output)
                nc_output = self.softmax(nc_output)
                ng_output = self.softmax(ng_output)
                nt_output = self.softmax(nt_output)

            return na_output, nc_output, ng_output, nt_output


class BasicConv2D(nn.Module):
    def __init__(self, input_channel, output_channel, kernel_size=3, strides=1, padding=1):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv2d(input_channel, output_channel, kernel_size=kernel_size, stride=strides, padding=padding,
                      bias=False),
            nn.BatchNorm2d(output_channel),
            nn.ReLU(inplace=True))

    def forward(self, x):
        return self.conv(x)


class ResNet(nn.Module):

    def __init__(self, block=BasicBlock, num_block=1, platform='ont', num_classes=3):
        super().__init__()

        in_channels = param.channel_size
        conv1_channel_size = 64
        conv2_channel_size = 96
        conv3_channel_size = 128
        fc1_layer_size = 160
        fc2_layer_size = 128
        dropout_rate = 0.3

        self.conv1 = BasicConv2D(input_channel=in_channels, output_channel=conv1_channel_size, strides=2)

        depth_scale_size = cal_scale(input_size=param.matrix_depth_dict[platform], layers=3)
        width_scale_size = cal_scale(input_size=param.no_of_positions, layers=3)

        self.conv1_x = self._make_layer(block, conv1_channel_size)

        self.conv2 = BasicConv2D(input_channel=conv1_channel_size, output_channel=conv2_channel_size, strides=2)
        self.conv2_x = self._make_layer(block, conv2_channel_size)

        self.conv3 = BasicConv2D(input_channel=conv2_channel_size, output_channel=conv3_channel_size, strides=2)
        self.conv3_x = self._make_layer(block, conv3_channel_size)

        self.dropout = nn.Dropout(p=dropout_rate)
        self.selu = nn.SELU()

        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(conv3_channel_size * depth_scale_size * width_scale_size, fc1_layer_size)
        self.fc2 = nn.Linear(fc1_layer_size, fc2_layer_size)
        self.fc3 = nn.Linear(fc2_layer_size, num_classes)

    def _make_layer(self, block, out_channels, stride=1):
        strides = [stride]
        layers = []
        for stride in strides:
            layers.append(block(out_channels, out_channels, stride))
            self.in_channels = out_channels * block.expansion

        return nn.Sequential(*layers)

    def forward(self, x):
        output = self.conv1(x)
        output = self.conv1_x(output)
        output = self.conv2(output)
        output = self.conv2_x(output)
        output = self.conv3(output)
        output = self.conv3_x(output)

        output = self.flatten(output)
        output = output.view(output.size(0), -1)
        output = self.dropout(output)
        output = self.selu(self.dropout(self.fc1(output)))
        output = self.selu(self.dropout(self.fc2(output)))
        output = self.selu(self.fc3(output))

        return output
