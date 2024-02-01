<div align="center">
    <img src="images/ClairS-TO_icon.png" width="200" alt="ClairS-TO">
</div>

# ClairS-TO - a deep-learning method for tumor-only somatic SNV calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng, Lei Chen  
Email: {rbluo,zxzheng,lchen}@cs.hku.hk 

------

## Introduction

ClairS-TO (Somatic Tumor-Only) is a tool in the Clair series to support long-read somatic single nucleotide variation (SNV) calling with only tumor samples available.

Without a normal sample, non-somatic noises cannot be identified by finding common signals between a paired tumor and normal. The variant caller itself needs to be more proficient in telling noises from somatic signals.

In ClairS-TO, we use an ensemble of two neural networks with opposite objectives. With the same input, an Affirmative NN determines how likely a candidate is a somatic variant - P(*Y<sub>Aff</sub>*), and a Negational NN determines how likely a candidate is NOT a somatic variant - P(*Y<sub>Neg</sub>*). A conditional probability P(*Y<sub>Aff</sub>* | *Y<sub>Neg</sub>*) that determines how likely a candidate is a somatic variant given the probability that the candidate is not a somatic variant is calculated from the probability of both networks. A somatic variant candidate that doesn't look like a noise usually has a high P(*Y<sub>Aff</sub>*) but a low P(*Y<sub>Neg</sub>*), while a somatic variant candidate that can also be a noise can have both a high P(*Y<sub>Aff</sub>*) and a high P(*Y<sub>Neg</sub>*).

Below is a workflow of ClairS-TO.
![ClairS-TO Workflow](images/ClairS-TO_architecture.png)

Like other tumor-only somatic variant callers, ClairS-TO accepts panel of normals (PoNs) as input to remove non-somatic variants.

For somatic variant calling using paired tumor/normal samples, please try [ClairS](https://github.com/HKU-BAL/ClairS).

------

## Contents
- [Latest Updates](#latest-updates)
- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Singularity](#option-2-singularity)
  - [Option 3. Build a micromamba (or anaconda) virtual environment](#option-3-build-a-micromamba-or-anaconda-virtual-environment)
  - [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
- [Quick Demo](#quick-demo)
- [Pre-trained Models](#pre-trained-models)
- [Usage](#usage)
- [Tagging non-somatic variant using panel of normals](#tagging-non-somatic-variant-using-panel-of-normals)
- [Disclaimer](#disclaimer)

------

## Latest Updates

*v0.0.2 (Jan. 26, 2024)*: 1. Added ONT Guppy 5kHz HAC (`-p ont_r10_guppy_hac_5khz`) and Dorado 4kHz HAC (`-p ont_r10_dorado_hac_4khz`) models, check [here](#pre-trained-models) for more details. 2. Added `FAU`, `FCU`, `FGU`, `FTU`, `RAU`, `RCU`, `RGU`, and `RTU` tags for the count of forward/reverse strand reads supporting A/C/G/T. 3. Revamped the way how panel of normals (PoNs) are inputted. Population databases are also considered as PoNs, and users can disable default population databases and add multiple other PoNs. 4. Added `file` and `md5` information of the PoNs to the VCF output header. 5. Enabled somatic variant calling in sex chromosomes. 6. Fixed an issue that misses PoNs tagging for low-quality variants.

*v0.0.1 (Dec. 4, 2023)*: Initial release for early access.

---

## Quick Demo

- Oxford Nanopore (ONT) [Q20+](https://nanoporetech.com/q20plus-chemistry) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- PacBio HiFi Revio data as input, see [PacBio HiFi Quick Demo](docs/pacbio_hifi_quick_demo.md).
- Illumina NGS data as input, see [Illumina Quick Demo](docs/illumina_quick_demo.md).

### Quick start

After following [installation](#installation), you can run ClairS-TO with one command:

```bash
./run_clairs_to -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy_sup_4khz

## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

------

## Pre-trained Models

ClairS-TO trained both Affirmative and Negational models using GIAB samples, and carry on benchmarking on HCC1395 tumor sample dataset. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|  Platform   |        Model name         |      Chemistry /Instruments      | Basecaller | Latest update | Option (`-p/--platform`)  |   Reference   |  Aligner   |
|:-----------:|:-------------------------:|:--------------------------------:|:----------:|:-------------:|:-------------------------:|:-------------:|:----------:|
|     ONT     | r1041_e82_400bps_sup_v420 |          R10.4.1, 5khz           | Dorado SUP | Nov. 10, 2023 | `ont_r10_dorado_sup_5khz` | GRCh38_no_alt |  Minimap2  |
|     ONT     | r1041_e82_400bps_sup_v410 |          R10.4.1, 4khz           | Dorado SUP | Nov. 10, 2023 | `ont_r10_dorado_sup_4khz` | GRCh38_no_alt |  Minimap2  |
|     ONT     | r1041_e82_400bps_hac_v410 |          R10.4.1, 4khz           | Dorado HAC | Jan. 19, 2024 | `ont_r10_dorado_hac_4khz` | GRCh38_no_alt |  Minimap2  |
|     ONT     | r1041_e82_400bps_sup_g615 |          R10.4.1, 4khz           | Guppy6 SUP | Nov. 10, 2023 | `ont_r10_guppy_sup_4khz`  | GRCh38_no_alt |  Minimap2  |
|     ONT     | r1041_e82_400bps_hac_g657 |          R10.4.1, 5khz           | Guppy6 HAC | Jan. 21, 2024 | `ont_r10_guppy_hac_5khz`  | GRCh38_no_alt |  Minimap2  |
|  Illumina   |           ilmn            |          NovaSeq/HiseqX          |     -      | Nov. 10, 2023 |          `ilmn`           |    GRCh38     |  BWA-MEM   |
| PacBio HiFi |        hifi_revio         | Revio with SMRTbell prep kit 3.0 |     -      | Nov. 10, 2023 |       `hifi_revio`        | GRCh38_no_alt |  Minimap2  |

------


## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/clairs-to). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs-to:latest \
  /opt/bin/run_clairs_to \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_hac_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy_sup_4khz, ont_r10_guppy_hac_5khz, ilmn, hifi_revio}
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
singularity pull docker://hkubal/clairs-to:latest

# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clairs-to_latest.sif \
  hkubal/clairs-to:latest \
  /opt/bin/run_clairs_to \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \      ## use your tumor bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_hac_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy_sup_4khz, ont_r10_guppy_hac_5khz, ilmn, hifi_revio}
  --output_dir ${OUTPUT_DIR} \                 ## output path prefix
  --conda_prefix /opt/micromamba/envs/clairs-to
```

### Option 3. Build a micromamba (or anaconda) virtual environment

Check here to install the tools step by step.

**Use micromamba (recommended)**:

Please install micromamba using the official [guide](https://mamba.readthedocs.io/en/latest/micromamba-installation.html) or using the commands below:

```bash
wget -O linux-64_micromamba-1.5.1-2.tar.bz2 https://micro.mamba.pm/api/micromamba/linux-64/latest
mkdir micromamba
tar -xvjf linux-64_micromamba-1.5.1-2.tar.bz2 -C micromamba
cd micromamba
./bin/micromamba shell init -s bash -p .
source ~/.bashrc
```

**Or use anaconda**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install ClairS-TO using micromamba step by step:**

```bash
# create and activate an environment named clairs-to
# install pypy and packages in the environment
# for micromamba
micromamba create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3-illumina bcftools einops python=3.9.0 -y
micromamba activate clairs-to

## for anaconda 
#conda create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3-illumina bcftools einops python=3.9.0 -y
#source activate clairs-to

git clone https://github.com/HKU-BAL/ClairS-TO.git
cd ClairS-TO

# make sure in clairs-to environment
# download pre-trained models
echo ${CONDA_PREFIX}
mkdir -p ${CONDA_PREFIX}/bin/clairs-to_models
mkdir -p ${CONDA_PREFIX}/bin/clairs-to_databases
wget http://www.bio8.cs.hku.hk/clairs-to/models/clairs-to_models.tar.gz
wget http://www.bio8.cs.hku.hk/clairs-to/databases/clairs-to_databases.tar.gz
tar -zxvf clairs-to_models.tar.gz -C ${CONDA_PREFIX}/bin/clairs-to_models/
tar -zxvf clairs-to_databases.tar.gz -C ${CONDA_PREFIX}/bin/clairs-to_databases/

./run_clairs_to --help
```

### Option 4. Docker Dockerfile

This is the same as Option 1 except that you are building a docker image yourself. Please refer to Option 1 for usage. 

```bash
git clone https://github.com/HKU-BAL/ClairS-TO.git
cd ClairS-TO

# build a docker image named hkubal/clairs-to:latest
# might require docker authentication to build docker image
docker build -f ./Dockerfile -t hkubal/clairs-to:latest .

# run the docker image like Option 1
docker run -it hkubal/clairs-to:latest /opt/bin/run_clairs_to --help
```

------

## Usage

### General Usage

```bash
./run_clairs_to \
  --tumor_bam_fn ${INPUT_DIR}/tumor.bam \    ## use your tumor bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \             ## use your reference file name here
  --threads ${THREADS} \                     ## maximum threads to be used
  --platform ${PLATFORM} \                   ## options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_hac_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy_sup_4khz, ont_r10_guppy_hac_5khz, ilmn, hifi_revio}
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
  -p, --platform PLATFORM           Select the sequencing platform of the input. Possible options {ont_r10_dorado_sup_4khz, ont_r10_dorado_hac_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy_sup_4khz, ont_r10_guppy_hac_5khz, ilmn, hifi_revio}.
```

**Miscellaneous parameters:**

```bash
  --pileup_affirmative_model_path PILEUP_AFFIRMATIVE_MODEL_PATH                                                                                                                                                
                        Specify the path to your own tumor-only somatic calling pileup affirmative model.                                                                                                                                                                                     
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
                        Minimal SNV AF required for a variant to be called. Decrease SNV_MIN_AF might increase a bit of sensitivity, but in trade of precision, speed, and accuracy. Default: 0.05.
  --min_coverage MIN_COVERAGE
                        Minimal coverage required for a variant to be called. Default: 4.
  --panel_of_normals PANEL_OF_NORMALS                                                                
                        The path of the panel of normals (PoNs) used for tagging non-somatic variants. Split by ',' for multiple PoNs. Default: if not specified, provided 'gnomad.r2.1.af-ge-0.001.sites.vcf.gz', 'dbsnp.b138.non-somatic.sites.vcf.gz', and '1000g-pon.sites.vcf.gz' PoNs will be included.
  --disable_nonsomatic_tagging
                        Disable non-somatic variants tagging. Default: enable non-somatic variants tagging.
  --chunk_size CHUNK_SIZE                           
                        The size of each chuck for parallel processing. Default: 5000000.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file. Default: SAMPLE.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output VCF filename. Default: output.
  --remove_intermediate_dir
                        Remove the intermediate directory before finishing to save disk space.
  --include_all_ctgs    Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}.
  --print_ref_calls     Show reference calls (0/0) in VCF file.
  --disable_print_nonsomatic_calls
                        Disable printing non-somatic calls. Default: enable non-somatic calls printing.
  -d, --dry_run         Print the commands that will be run.
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
              
```

#### Call SNVs in one or multiple chromosomes using the `-C/--ctg_name` parameter

```bash
./run_clairs_to -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy_sup_4khz -C chr21,chr22
```

#### Call SNVs in one specific region using the `-r/--region` parameter

```bash
./run_clairs_to -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy_sup_4khz -r chr20:1000000-2000000
```

#### Call SNVs at interested variant sites (genotyping) using the `-G/--genotyping_mode_vcf_fn` parameter

```bash
./run_clairs_to -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy_sup_4khz -G input.vcf
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
./run_clairs_to -T tumor.bam -R ref.fa -o output -t 8 -p ont_r10_guppy_sup_4khz -B input.bed
```

------

## Tagging non-somatic variant using panel of normals

ClairS-TO by default tags variants if they exist in provided panel of normals (PoNs, i.e., `gnomad.r2.1.af-ge-0.001.sites.vcf.gz`, `dbsnp.b138.non-somatic.sites.vcf.gz`, and `1000g-pon.sites.vcf.gz`), and pass the filters listed in the table below. 

Users can also use their own PoNs for tagging using the `--panel_of_normals` option. 

Particularly, if the `--panel_of_normals` option is not specified, the three default PoNs will be included. And if users want to use all/part/none of the default PoNs as well as their own PoNs, corresponding file paths of the default PoNs (i.e., `${CONDA_PREFIX}/bin/clairs-to_databases/gnomad.r2.1.af-ge-0.001.sites.vcf.gz`, `${CONDA_PREFIX}/bin/clairs-to_databases/dbsnp.b138.non-somatic.sites.vcf.gz`, and `${CONDA_PREFIX}/bin/clairs-to_databases/1000g-pon.sites.vcf.gz`), and their own PoNs, should be included in the `--panel_of_normals` option, split by `,`. 

| Default PoNs |                                        URL                                         | Source |                                                    Source URL                                                    |      Last visited        | Total #Variants  |        Filters        | #Variants used for tagging | Remaining Columns in the input   |
|:-------------:|:----------------------------------------------------------------------------------:|:------:|:----------------------------------------------------------------------------------------------------------------:|:------------------------:|:----------------:|:---------------------:|:--------------------------:|:--------------------------------:|
|   PoN 1    | http://www.bio8.cs.hku.hk/clairs-to/databases/gnomad.r2.1.af-ge-0.001.sites.vcf.gz |  GATK gnomAD  |            https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz            | July 10, 2023 PM10∶34∶07 |   268,225,276    | Sites with AF ≥ 0.001 |         35,551,905         |     #CHROM  POS ID  REF ALT      |
|   PoN 2   | http://www.bio8.cs.hku.hk/clairs-to/databases/dbsnp.b138.non-somatic.sites.vcf.gz  |  GATK dbSNP  | https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf | July 10, 2023 PM10∶42∶22 |    60,691,395    |   Non-Somatic sites   |         60,683,019         |     #CHROM  POS ID  REF ALT      |
|   PoN 3   |        http://www.bio8.cs.hku.hk/clairs-to/databases/1000g-pon.sites.vcf.gz        |  GATK 1000G PoN  |              https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz               | July 10, 2023 PM10∶31∶32 |    2,609,566     |       All sites       |         2,609,566          |     #CHROM  POS ID  REF ALT      |

------

## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.
