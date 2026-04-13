#!/bin/bash

# Load conda + activate environment
eval "$(/home/chese-u/miniconda/bin/conda shell.bash hook)"
conda activate bio
# Script to filter and summarize variants (WSL version)

# directories - להתאים לנתיבים שלך ב־demo2
ref="$HOME/demo2/supporting_files/hg38/hg38.fa"
results="$HOME/demo2/results"

############################################
# חלק א: סינון ויצירת analysis-ready VCFs #
############################################

# אפשר להשאיר את הבלוק הזה ב־if false אם כבר הרצת ולא רוצה לדרוס קבצים.
# אם תרצה להריץ שוב, שנה ל: if true

if false
then

  # Filter SNPs
  gatk VariantFiltration \
    -R "${ref}" \
    -V "${results}/raw_snps.vcf" \
    -O "${results}/filtered_snps.vcf" \
    -filter-name "QD_filter"         -filter "QD < 2.0" \
    -filter-name "FS_filter"         -filter "FS > 60.0" \
    -filter-name "MQ_filter"         -filter "MQ < 40.0" \
    -filter-name "SOR_filter"        -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter"  -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

  # Filter INDELs
  gatk VariantFiltration \
    -R "${ref}" \
    -V "${results}/raw_indels.vcf" \
    -O "${results}/filtered_indels.vcf" \
    -filter-name "QD_filter"   -filter "QD < 2.0" \
    -filter-name "FS_filter"   -filter "FS > 200.0" \
    -filter-name "SOR_filter"  -filter "SOR > 10.0" \
    -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

  # Select only variants that pass FILTER (PASS)
  gatk SelectVariants \
    --exclude-filtered \
    -V "${results}/filtered_snps.vcf" \
    -O "${results}/analysis-ready-snps.vcf"

  gatk SelectVariants \
    --exclude-filtered \
    -V "${results}/filtered_indels.vcf" \
    -O "${results}/analysis-ready-indels.vcf"

  # להוציא גם גנוטיפים שנכשלו בפילטרי DP/GQ
  grep -v -E "DP_filter|GQ_filter" "${results}/analysis-ready-snps.vcf" \
    > "${results}/analysis-ready-snps-filteredGT.vcf"

  grep -v -E "DP_filter|GQ_filter" "${results}/analysis-ready-indels.vcf" \
    > "${results}/analysis-ready-indels-filteredGT.vcf"

fi

#####################################
# חלק ב: יצוא לטבלת SNPs (ללא R)   #
#####################################

gatk VariantsToTable \
  -V "${results}/analysis-ready-snps-filteredGT.vcf" \
  -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
  -F AC -F AN -F DP -F AF \
  -O "${results}/output_snps.table"

