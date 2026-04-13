#!/bin/bash
set -euo pipefail

#####################################
# 0. activate conda env "bio"
#####################################

# load conda into this shell
eval "$(/home/chese-u/miniconda/bin/conda shell.bash hook)"
conda activate bio

echo "=== Germline variant calling demo on WSL ==="
date
echo

#####################################
# 1. directories and paths
#####################################

BASE="/home/chese-u/demo2"

reads="$BASE/reads"
aligned="$BASE/aligned_reads"
results="$BASE/results"
data="$BASE/data"
ref_dir="$BASE/supporting_files/hg38"
log="$BASE/pipeline.log"

ref="$ref_dir/hg38.fa"
known_sites="$ref_dir/Homo_sapiens_assembly38.dbsnp138.vcf"

mkdir -p "$reads" "$aligned" "$results" "$data" "$ref_dir"

echo "Using base dir: $BASE"
echo "Log file: $log"
echo

#####################################
# helper function
#####################################

step() {
  echo
  echo "=== $1 ==="
  echo "=== $1 ===" >> "$log"
  date | tee -a "$log"
}

#####################################
# STEP 0: reference and known sites
#####################################

step "STEP 0: prepare reference and known sites"

# download hg38 reference if missing
if [ ! -f "$ref" ]; then
  echo "Downloading hg38 reference..." | tee -a "$log"
  wget -O "$ref_dir/hg38.fa.gz" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip "$ref_dir/hg38.fa.gz"
else
  echo "Reference fasta already exists: $ref" | tee -a "$log"
fi

# samtools index .fai
if [ ! -f "$ref.fai" ]; then
  echo "Building fasta index (.fai) with samtools faidx..." | tee -a "$log"
  samtools faidx "$ref"
else
  echo "Fasta index already exists: $ref.fai" | tee -a "$log"
fi

# GATK sequence dictionary .dict
if [ ! -f "$ref_dir/hg38.dict" ]; then
  echo "Creating sequence dictionary with GATK..." | tee -a "$log"
  gatk CreateSequenceDictionary \
    -R "$ref" \
    -O "$ref_dir/hg38.dict"
else
  echo "Sequence dictionary already exists: $ref_dir/hg38.dict" | tee -a "$log"
fi

# download known sites for BQSR
if [ ! -f "$known_sites" ]; then
  echo "Downloading dbSNP known sites..." | tee -a "$log"
  wget -P "$ref_dir" \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
  wget -P "$ref_dir" \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
else
  echo "Known sites VCF already exists: $known_sites" | tee -a "$log"
fi

#####################################
# STEP 1: download reads + FastQC
#####################################

step "STEP 1: download reads and run FastQC"

read1="$reads/SRR062634_1.filt.fastq.gz"
read2="$reads/SRR062634_2.filt.fastq.gz"

if [ ! -f "$read1" ]; then
  echo "Downloading read 1..." | tee -a "$log"
  wget -P "$reads" \
    ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
else
  echo "Read 1 already exists: $read1" | tee -a "$log"
fi

if [ ! -f "$read2" ]; then
  echo "Downloading read 2..." | tee -a "$log"
  wget -P "$reads" \
    ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
else
  echo "Read 2 already exists: $read2" | tee -a "$log"
fi

# FastQC (לא מזיק להריץ שוב, אבל אפשר לבדוק HTML פעם אחת בלבד)
fastqc "$read1" -o "$reads"
fastqc "$read2" -o "$reads"

#####################################
# STEP 2: mapping with BWA MEM
#####################################

step "STEP 2: map to reference with BWA MEM"

# build BWA index once (.bwt זה אינדיקציה טובה)
if [ ! -f "$ref.bwt" ]; then
  echo "Building BWA index..." | tee -a "$log"
  bwa index "$ref"
else
  echo "BWA index already exists (found $ref.bwt)" | tee -a "$log"
fi

sam_out="$aligned/SRR062634.paired.sam"

if [ ! -f "$sam_out" ]; then
  echo "Running BWA MEM..." | tee -a "$log"
  bwa mem -t 4 \
    -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" \
    "$ref" \
    "$read1" \
    "$read2" \
    > "$sam_out"
else
  echo "SAM file already exists: $sam_out" | tee -a "$log"
fi

############################################
# STEP 3: MarkDuplicates + sort (GATK4)
############################################

step "STEP 3: mark duplicates and sort"

dedup_bam="$aligned/SRR062634.sorted.dedup.bam"
dedup_metrics="$aligned/SRR062634.dedup.metrics.txt"

if [ ! -f "$dedup_bam" ]; then
  echo "Running GATK MarkDuplicatesSpark..." | tee -a "$log"
  gatk MarkDuplicatesSpark \
    -I "$sam_out" \
    -O "$dedup_bam" \
    -M "$dedup_metrics"
else
  echo "Dedup BAM already exists: $dedup_bam" | tee -a "$log"
fi

#####################################
# STEP 4: BQSR
#####################################

step "STEP 4: base quality score recalibration (BQSR)"

recal_table="$data/recal_data.table"
bqsr_bam="$aligned/SRR062634.sorted.dedup.bqsr.bam"

if [ ! -f "$recal_table" ]; then
  echo "Running BaseRecalibrator..." | tee -a "$log"
  gatk BaseRecalibrator \
    -I "$dedup_bam" \
    -R "$ref" \
    --known-sites "$known_sites" \
    -O "$recal_table"
else
  echo "Recal table already exists: $recal_table" | tee -a "$log"
fi

if [ ! -f "$bqsr_bam" ]; then
  echo "Running ApplyBQSR..." | tee -a "$log"
  gatk ApplyBQSR \
    -I "$dedup_bam" \
    -R "$ref" \
    --bqsr-recal-file "$recal_table" \
    -O "$bqsr_bam"
else
  echo "BQSR BAM already exists: $bqsr_bam" | tee -a "$log"
fi

#####################################
# STEP 5: metrics and histogram
#####################################

step "STEP 5: collect alignment and insert size metrics"

align_metrics="$aligned/alignment_metrics.txt"
insert_metrics="$aligned/insert_size_metrics.txt"
insert_hist="$aligned/insert_size_histogram.pdf"

if [ ! -f "$align_metrics" ]; then
  gatk CollectAlignmentSummaryMetrics \
    R="$ref" \
    I="$bqsr_bam" \
    O="$align_metrics"
else
  echo "Alignment metrics already exist: $align_metrics" | tee -a "$log"
fi

if [ ! -f "$insert_metrics" ]; then
  gatk CollectInsertSizeMetrics \
    INPUT="$bqsr_bam" \
    OUTPUT="$insert_metrics" \
    HISTOGRAM_FILE="$insert_hist"
else
  echo "Insert size metrics already exist: $insert_metrics" | tee -a "$log"
fi

#####################################
# STEP 6: variant calling
#####################################

step "STEP 6: call variants with HaplotypeCaller"

raw_vcf="$results/raw_variants.vcf"
snps_vcf="$results/raw_snps.vcf"
indels_vcf="$results/raw_indels.vcf"

if [ ! -f "$raw_vcf" ]; then
  gatk HaplotypeCaller \
    -R "$ref" \
    -I "$bqsr_bam" \
    -O "$raw_vcf"
else
  echo "Raw variants VCF already exists: $raw_vcf" | tee -a "$log"
fi

if [ ! -f "$snps_vcf" ]; then
  gatk SelectVariants \
    -R "$ref" \
    -V "$raw_vcf" \
    --select-type SNP \
    -O "$snps_vcf"
else
  echo "SNP VCF already exists: $snps_vcf" | tee -a "$log"
fi

if [ ! -f "$indels_vcf" ]; then
  gatk SelectVariants \
    -R "$ref" \
    -V "$raw_vcf" \
    --select-type INDEL \
    -O "$indels_vcf"
else
  echo "INDEL VCF already exists: $indels_vcf" | tee -a "$log"
fi

echo
echo "Pipeline finished successfully."
date | tee -a "$log"
