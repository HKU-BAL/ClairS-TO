## PacBio HiFi Tumor-only Somatic Variant Calling Quick Demo
Here is a quick demo for the PacBio tumor-only somatic variant calling using HCC1395 tumor sample chromosome 17 data. The data was sequenced using PacBio latest Revio system.

```bash
Platform:          PacBio HiFi Revio
Sample:     	   HCC1395 tumor sample
Tumor coverage:    ~60x
Reference:         GRCh38_no_alt
Aligner:           pbmm2
Region:            chr17:80000000-80100000
```

**Download data**

```bash
# Parameters
INPUT_DIR="${HOME}/pacbio_hifi_quick_demo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/GRCh38_no_alt_chr17.fa
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/GRCh38_no_alt_chr17.fa.fai
# Tumor BAM
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/HCC1395_tumor_chr17_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/HCC1395_tumor_chr17_demo.bam.bai

# SEQC2 Truth VCF and BED
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz.tbi
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/pacbio_hifi/SEQC2_High-Confidence_Regions_v1.2_chr17.bed

REF="GRCh38_no_alt_chr17.fa"
TUMOR_BAM="HCC1395_tumor_chr17_demo.bam"
BASELINE_VCF_FILE_PATH="SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz"
BASELINE_BED_FILE_PATH="SEQC2_High-Confidence_Regions_v1.2_chr17.bed"
OUTPUT_VCF_FILE_PATH="output.vcf.gz"

```

#### Tumor-only somatic variant calling using docker pre-built image

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairsto:latest \
  /opt/bin/run_clairsto \
  --tumor_bam_fn ${INPUT_DIR}/${TUMOR_BAM} \
  --ref_fn ${INPUT_DIR}/${REF} \
  --threads 4 \
  --platform hifi_revio \
  --output_dir ${OUTPUT_DIR} \
  --region chr17:80000000-80100000
```

**Run [compare_vcf.py](src/compare.vcf) for benchmarking (optional)**

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairsto:latest \
  python3 /opt/bin/clairsto.py compare_vcf \
     --truth_vcf_fn ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
     --input_vcf_fn ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
     --bed_fn ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
     --output_dir ${OUTPUT_DIR}/benchmark \
     --input_filter_tag 'PASS' \
     --ctg_name chr17 \
     --ctg_start 80000000 \
     --ctg_end 80100000
```

**Expected output:**

|  Type   | Precision | Recall | F1-score | TP | FP | FN |
| :-----: |:---------:|:------:|:--------:|:--:|:--:|:--:|
| **SNV** |  0.9333   | 0.9655 |  0.9492  | 28 | 2  | 1  |


 **Or run [som.py]() for benchmarking (optional)**

```bash
# Need to restrict target BED regions for benchmarking
echo -e "chr17\t80000000\t80100000" > ${INPUT_DIR}/quick_demo.bed

docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/som.py \
    ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
    ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
    -T ${INPUT_DIR}/quick_demo.bed \
    -f ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
    -r ${INPUT_DIR}/${REF} \
    -o "${OUTPUT_DIR}/som" \
    -l chr17
```

**Run all commands above:**

```bash
cd ${HOME}
wget "https://raw.githubusercontent.com/HKU-BAL/clairsto/main/demo/pacbio_hifi_quick_demo.sh"
chmod +x pacbio_hifi_quick_demo.sh
./pacbio_hifi_quick_demo.sh
```

Check the results using `less ${HOME}/pacbio_hifi_quick_demo/output/output.vcf.gz`.

