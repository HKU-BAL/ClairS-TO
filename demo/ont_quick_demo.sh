# Parameters
INPUT_DIR="${HOME}/ont_quick_demo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/GRCh38_no_alt_chr17.fa
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/GRCh38_no_alt_chr17.fa.fai
# Tumor BAM
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/HCC1395_tumor_chr17_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/HCC1395_tumor_chr17_demo.bam.bai

# SEQC2 Truth VCF and BED
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz.tbi
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clairs/quick_demo/ont/SEQC2_High-Confidence_Regions_v1.2_chr17.bed

REF="GRCh38_no_alt_chr17.fa"
TUMOR_BAM="HCC1395_tumor_chr17_demo.bam"
BASELINE_VCF_FILE_PATH="SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz"
BASELINE_BED_FILE_PATH="SEQC2_High-Confidence_Regions_v1.2_chr17.bed"
OUTPUT_SNV_VCF_FILE_PATH="snv.vcf.gz"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs-to:latest \
  /opt/bin/run_clairs_to \
  --tumor_bam_fn ${INPUT_DIR}/${TUMOR_BAM} \
  --ref_fn ${INPUT_DIR}/${REF} \
  --threads 4 \
  --platform ont_r10_guppy_sup_4khz \
  --output_dir ${OUTPUT_DIR} \
  --region chr17:80000000-80100000

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clairs-to:latest \
  python3 /opt/bin/clairs_to.py compare_vcf \
     --truth_vcf_fn ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
     --input_vcf_fn ${OUTPUT_DIR}/${OUTPUT_SNV_VCF_FILE_PATH} \
     --bed_fn ${INPUT_DIR}/${BASELINE_BED_FILE_PATH} \
     --output_dir ${OUTPUT_DIR}/benchmark \
     --input_filter_tag 'PASS' \
     --ctg_name chr17 \
     --ctg_start 80000000 \
     --ctg_end 80100000
