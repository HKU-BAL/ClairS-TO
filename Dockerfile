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

# Example command:
# $ git clone https://github.com/HKU-BAL/ClairS-TO.git
# $ cd ClairS-TO
# $ docker build -f ./Dockerfile -t hkubal/clairs-to:latest .
# $ docker run -it hkubal/clairs-to:latest /opt/bin/run_clairs_to --help

FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/micromamba/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        g++ \
        time \
        libboost-graph-dev && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/bin

# install micromamba
RUN wget --quiet -O linux-64_micromamba-1.5.1-2.tar.bz2 https://micro.mamba.pm/api/micromamba/linux-64/latest && \
    tar -xvjf linux-64_micromamba-1.5.1-2.tar.bz2 && \
    mkdir -p /opt/micromamba/bin && \
    mv ./bin/micromamba /opt/micromamba/bin/micromamba && \
    /opt/micromamba/bin/micromamba shell init -s bash -p /opt/micromamba && \
    export MAMBA_ROOT_PREFIX=/opt/micromamba && \
    rm linux-64_micromamba-1.5.1-2.tar.bz2 && \
    rm -r info/ && \
    rm -r bin/ && \
    micromamba create -n clairs-to -c pytorch -c conda-forge -c bioconda clair3 bcftools einops tqdm pytorch torchinfo -y && \
    rm -rf /opt/micromamba/pkgs/* && \
    rm -rf /root/.cache/pip

ENV PATH /opt/micromamba/envs/clairs-to/bin:$PATH
ENV CONDA_DEFAULT_ENV clairs-to

RUN apt install curl zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev -y && \
    /opt/micromamba/envs/clairs-to/bin/python3 -m pip install scipy scikit-learn && \
    rm -rf /var/lib/apt/lists/*

COPY . .

RUN cd /opt/bin/src/realign && \
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    cd /opt/bin/src/verdict/allele_counter && chmod +x setup.sh && /bin/bash setup.sh /opt/bin/src/verdict/allele_counter && \
    wget http://www.bio8.cs.hku.hk/clairs-to/models/clairs-to_models.tar.gz	-P /opt/models && \
    wget http://www.bio8.cs.hku.hk/clairs-to/databases/clairs-to_databases.tar.gz -P /opt/databases && \
    wget http://www.bio8.cs.hku.hk/clairs-to/cna_data/reference_files.tar.gz -P /opt/cna_data && \
    mkdir -p /opt/micromamba/envs/clairs-to/bin/clairs-to_models && \
    mkdir -p /opt/micromamba/envs/clairs-to/bin/clairs-to_databases && \
    mkdir -p /opt/micromamba/envs/clairs-to/bin/clairs-to_cna_data && \
    tar -zxvf /opt/models/clairs-to_models.tar.gz -C /opt/micromamba/envs/clairs-to/bin/clairs-to_models && \
    tar -zxvf /opt/databases/clairs-to_databases.tar.gz -C /opt/micromamba/envs/clairs-to/bin/clairs-to_databases && \
    tar -zxvf /opt/cna_data/reference_files.tar.gz -C /opt/micromamba/envs/clairs-to/bin/clairs-to_cna_data && \
    rm /opt/models/clairs-to_models.tar.gz && \
    rm /opt/databases/clairs-to_databases.tar.gz && \
    rm -rf /opt/cna_data/reference_files.tar.gz && \
    echo 'will cite' | parallel --citation || true \
    echo "micromamba activate clairs-to" > ~/.bashrc
