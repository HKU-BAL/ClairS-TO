<!-- # ClairS-TO - a deep-learning method for tumor-only long-read somatic small variant calling -->
# ClairS-TO

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng, Lei Chen  
Email: rbluo@cs.hku.hk, zxzheng@cs.hku.hk, lchen@cs.hku.hk 

------

## Introduction

ClairS-TO is a tumor-only somatic variant caller without a paired normal sample and primarily for ONT long-read.
It calculates the probability of a somatic variant candidate using **Bayes' Theorem**, **Affirmational Neural Network** and **Negational Neural Network**.  

Specifically, the Affirmational Neural Network predicts the somatic variant candidate being an A, C, G or T, while the Negational Neural Network predicts it **NOT** being an A, C, G or T. 
In addition, the Bayes' Theorem is adopted to combine the predictions from the two Networks for final output.

Particularly, genetic databases (e.g., gnomAD, dbSNP, and 1000G PoN) are utilized to tag germline variants.

------

## Contents
- [Latest Updates](#latest-updates)
- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Singularity](#option-2-singularity)
  - [Option 3. Build an anaconda virtual environment](#option-3-build-an-anaconda-virtual-environment)
  - [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
- [Quick Demo](#quick-demo)
- [Pre-trained Models](#pre-trained-models)
- [Usage](#usage)
- [Disclaimer](#disclaimer)

------

## Latest Updates

*v0.0.1 (Nov., 2023)*: Initial release for early access.

---

## Quick Demo

- Oxford Nanopore (ONT) [Q20+](https://nanoporetech.com/q20plus-chemistry) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- PacBio HiFi [Revio](https://www.pacb.com/revio/) data as input, see [PacBio HiFi Quick Demo](docs/pacbio_hifi_quick_demo.md).
- Illumina NGS data as input, see [Illumina Quick Demo](docs/illumina_quick_demo.md).

### Quick start

After following [installation](#installation), you can run ClairS-TO with one command:

```bash
./run_clairsto -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy

## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

------

## Pre-trained Models

ClairS-TO trained both Affirmational and Negational models using GIAB samples, and carry on benchmarking on HCC1395 tumor sample dataset. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|  Platform   |        Model name         |      Chemistry /Instruments      | Basecaller | Option (`-p/--platform`) |   Reference   | Aligner  |
| :---------: |:-------------------------:|:--------------------------------:|:----------:| :-----------: | :------: | ----------- |
| ONT | r1041_e82_400bps_sup_v420 |          R10.4.1, 5khz           |   Dorado   | `ont_r10_dorado_5khz` | GRCh38_no_alt | Minimap2 |
| ONT | r1041_e82_400bps_sup_v410 |          R10.4.1, 4khz           |   Dorado   | `ont_r10_dorado_4khz` | GRCh38_no_alt | Minimap2 |
| ONT | r1041_e82_400bps_sup_g615 |              R10.4.1             |   Guppy6   |   `ont_r10_guppy` | GRCh38_no_alt | Minimap2 |
|  Illumina   |           ilmn            |          NovaSeq/HiseqX          |     -      |          `ilmn`          |    GRCh38     | BWA-MEM  |
| PacBio HIFI |        hifi_revio         | Revio with SMRTbell prep kit 3.0 |     -      | `hifi_revio` | GRCh38_no_alt | Minimap2 |

------


## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/clairsto). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairsto:latest \
  /opt/bin/run_clairsto \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_4khz, ont_r10_dorado_5khz, ont_r10_guppy, ilmn, hifi_revio}
  --output_dir ${OUTPUT_DIR}                   ## output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in singularity. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
mkdir -p ${OUTPUT_DIR}

conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clairsto:latest

# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clairsto_latest.sif \
  hkubal/clairsto:latest \
  /opt/bin/run_clairsto \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_4khz, ont_r10_dorado_5khz, ont_r10_guppy, ilmn, hifi_revio}
  --output_dir ${OUTPUT_DIR} \                 ## output path prefix
  --conda_prefix /opt/conda/envs/clairsto
```

### Option 3. Build an anaconda virtual environment

Check here to install the tools step by step.

**Anaconda install**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install ClairS-TO using anaconda step by step:**

```bash
# create and activate an environment named clairsto
# install pypy and packages in the environemnt
conda create -n clairsto -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3-illumina python=3.9.0 -y
source activate clairsto

git clone https://github.com/HKU-BAL/ClairS-TO.git
cd ClairS-TO

# make sure in conda environment
# download pre-trained models
echo ${CONDA_PREFIX}
mkdir -p ${CONDA_PREFIX}/bin/clairsto_models
wget http://www.bio8.cs.hku.hk/clairsto/models/clairsto_models.tar.gz
tar -zxvf clairsto_models.tar.gz -C ${CONDA_PREFIX}/bin/clairsto_models/

./run_clairsto --help
```

### Option 4. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
git clone https://github.com/HKU-BAL/ClairS-TO.git
cd ClairS-TO

# build a docker image named hkubal/clairsto:latest
# might require docker authentication to build docker image
docker build -f ./Dockerfile -t hkubal/clairsto:latest .

# run the docker image like option 1
docker run -it hkubal/clairsto:latest /opt/bin/run_clairsto --help
```

------

## Usage

### General Usage

```bash
./run_clairsto \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \    ## use your tumor bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \             ## use your reference file name here
  --threads ${THREADS} \                     ## maximum threads to be used
  --platform ${PLATFORM} \                   ## options: {ont_r10_dorado_4khz, ont_r10_dorado_5khz, ont_r10_guppy, ilmn, hifi_revio}
  --output_dir ${OUTPUT_DIR}                 ## output path prefix
 
## Final output file: ${OUTPUT_DIR}/output.vcf.gz
```

### Options

**Required parameters:**

```bash
  -T, --tumor_bam_fn TUMOR_BAM_FN   Tumor BAM file input. The input file must be samtools indexed.
  -R, --ref_fn FASTA                Reference file input. The input file must be samtools indexed.
  -o, --output_dir OUTPUT_DIR       VCF output directory.
  -t, --threads THREADS             Max threads to be used.
  -p, --platform PLATFORM           Select the sequencing platform of the input. Possible options {ont_r10_dorado_4khz, ont_r10_dorado_5khz, ont_r10_guppy, ilmn, hifi_revio}.
```

**Miscellaneous parameters:**

```bash
  --pileup_affirmational_model_path PILEUP_AFFIRMATIONAL_MODEL_PATH                                                                                                                                                
                        Specify the path to your own tumor-only somatic calling pileup affirmational model.                                                                                                                                                                                     
  --pileup_negational_model_path PILEUP_NEGATIONAL_MODEL_PATH
                        Specify the path to your own tumor-only somatic calling pileup negational model.
  -c CTG_NAME, --ctg_name CTG_NAME
                        The name of the contigs to be processed. Split by ',' for multiple contigs. Default: all contigs will be processed.
  -r REGION, --region REGION
                        A region to be processed. Format: `ctg_name:start-end` (start is 1-based).
  -b BED_FN, --bed_fn BED_FN
                        Path to a BED file. Call variants only in the provided BED regions.
  -G GENOTYPING_MODE_VCF_FN, --genotyping_mode_vcf_fn GENOTYPING_MODE_VCF_FN
                        VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided.
  -H HYBRID_MODE_VCF_FN, --hybrid_mode_vcf_fn HYBRID_MODE_VCF_FN
                        Enable hybrid calling mode that combines the de novo calling results and genotyping results at the positions in the VCF file given.
  -q QUAL, --qual QUAL  If set, variants with >QUAL will be marked as PASS, or LowQual otherwise.
  --snv_min_af SNV_MIN_AF                           
                        Minimal SNV AF required for a variant to be called. Decrease SNV_MIN_AF might increase a bit of sensitivity, but in trade of precision, speed and accuracy. Default: 0.05.
  --min_coverage MIN_COVERAGE
                        Minimal coverage required for a variant to be called. Default: 4.
  --chunk_size CHUNK_SIZE                           
                        The size of each chuck for parallel processing. Default: 5000000.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file. Default: SAMPLE.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output VCF filename. Default: output.
  --remove_intermediate_dir
                        Remove intermediate directory before finishing to save disk space.
  --include_all_ctgs    Call variants on all contigs, otherwise call in chr{1..22} and {1..22}.
  --print_ref_calls     Show reference calls (0/0) in VCF file.
  --disable_print_germline_calls
                        Disable print germline calls. Default: enable germline calls printing.
  -d, --dry_run         Print the commands that will be ran.
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
  --disable_germline_tagging
                        Disable germline tagging to tag germline calls. Default: enable germline tagging.
  --disable_gnomad_tagging
                        Disable gnomAD database resource to tag germline calls. Default: enable gnomad tagging.
  --disable_pon_tagging
                        Disable 1000G PoN database resource to tag germline calls. Default: enable pon tagging.
  --disable_dbsnp_tagging
                        Disable dbSNP database resource to tag germline calls. Default: enable dbsnp tagging.
  --use_own_pon_resource USE_OWN_PON_RESOURCE
                        Use own PoN resource to tag germline calls.                  
```

#### Call SNVs in one or mutiple chromosomes using the `-C/--ctg_name` parameter

```bash
./run_clairsto -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -C chr21,chr22
```

#### Call SNVs in one specific region using the `-r/--region` parameter

```bash
./run_clairsto -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -r chr20:1000000-2000000
```

#### Call SNVs at interested variant sites (genotyping) using the `-G/--genotyping_mode_vcf_fn` parameter

```bash
./run_clairsto -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -G input.vcf
```

#### Call SNVs in the BED regions using the `-B/--bed_fn` parameter

We highly recommended using BED file to define multiple regions of interest like:

```shell
echo -e "${CTG1}\t${START_POS_1}\t${END_POS_1}" > input.bed
echo -e "${CTG2}\t${START_POS_2}\t${END_POS_2}" >> input.bed
...
```

Then:

```bash
./run_clairsto -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy -B input.bed
```

------


## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.
