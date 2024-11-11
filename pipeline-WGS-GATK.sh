#!/bin/bash

# Full pipeline script for genetic variant calling, filtering, and annotation using GATK4
# This script demonstrates a complete bioinformatics pipeline for genomic analysis

# Define directories and file paths (adjust these as needed for your environment)
ref="reference_genome/hg38.fa"                         # Reference genome path (after download)
reads_dir="input_reads"                                # Directory for input reads
results="output_directory/results"                     # Directory to store results
funcotator_data="funcotator_dataSources"               # Funcotator data sources path

# Create directories if they don’t exist
mkdir -p ${results} ${reads_dir} ${funcotator_data} reference_genome

# --------------------------------------
# PART 0: Download Required Files
# --------------------------------------

# Download reference genome (adjust source if needed)
if [ ! -f ${ref} ]; then
    echo "Downloading reference genome..."
    wget -O ${ref} "https://genome-source-url/hg38.fa"
    # Download associated reference files like .fai and .dict if necessary
    wget -O reference_genome/hg38.fa.fai "https://genome-source-url/hg38.fa.fai"
    wget -O reference_genome/hg38.dict "https://genome-source-url/hg38.dict"
fi

# Download known sites for BQSR (e.g., dbSNP and known indels)
if [ ! -f "reference_genome/dbsnp.vcf" ]; then
    echo "Downloading dbSNP..."
    wget -O reference_genome/dbsnp.vcf "https://resource-url/dbsnp.vcf"
fi
if [ ! -f "reference_genome/known_indels.vcf" ]; then
    echo "Downloading known indels..."
    wget -O reference_genome/known_indels.vcf "https://resource-url/known_indels.vcf"
fi

# Download input reads for demonstration
if [ ! -f "${reads_dir}/sample_R1.fastq" ] || [ ! -f "${reads_dir}/sample_R2.fastq" ]; then
    echo "Downloading sample input reads..."
    wget -O ${reads_dir}/sample_R1.fastq "https://reads-source-url/sample_R1.fastq"
    wget -O ${reads_dir}/sample_R2.fastq "https://reads-source-url/sample_R2.fastq"
fi

# Download Funcotator data sources (adjust URL as needed)
if [ ! -d ${funcotator_data} ]; then
    echo "Downloading Funcotator data sources..."
    wget -O ${funcotator_data}/funcotator_data.tar.gz "https://funcotator-data-url/funcotator_dataSources.tar.gz"
    tar -xzvf ${funcotator_data}/funcotator_data.tar.gz -C ${funcotator_data}
fi

# --------------------------------------
# PART 1: Preprocessing and Variant Calling
# --------------------------------------

# Step 1: Align reads to the reference genome using BWA-MEM
bwa mem -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" ${ref} ${reads_dir}/sample_R1.fastq ${reads_dir}/sample_R2.fastq > ${results}/aligned_reads.sam

# Step 2: Convert SAM to BAM, sort, and index using samtools
samtools view -S -b ${results}/aligned_reads.sam > ${results}/aligned_reads.bam
samtools sort ${results}/aligned_reads.bam -o ${results}/sorted_reads.bam
samtools index ${results}/sorted_reads.bam

# Step 3: Mark duplicates using GATK
gatk MarkDuplicates -I ${results}/sorted_reads.bam -O ${results}/dedup_reads.bam -M ${results}/marked_dup_metrics.txt
samtools index ${results}/dedup_reads.bam

# Step 4: Base quality score recalibration (BQSR) using GATK
gatk BaseRecalibrator -R ${ref} -I ${results}/dedup_reads.bam --known-sites reference_genome/dbsnp.vcf --known-sites reference_genome/known_indels.vcf -O ${results}/recal_data.table
gatk ApplyBQSR -R ${ref} -I ${results}/dedup_reads.bam --bqsr-recal-file ${results}/recal_data.table -O ${results}/recal_reads.bam

# Step 5: Variant calling using GATK HaplotypeCaller
gatk HaplotypeCaller -R ${ref} -I ${results}/recal_reads.bam -O ${results}/raw_variants.vcf

# Split SNPs and INDELs into separate files for further filtering
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type-to-include SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type-to-include INDEL -O ${results}/raw_indels.vcf

# --------------------------------------
# PART 2: Filter and Annotate Variants
# --------------------------------------

# Part 1: Filter SNPs with specific quality criteria
gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/raw_snps.vcf \
    -O ${results}/filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

# Part 2: Filter INDELs with specific quality criteria
gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/raw_indels.vcf \
    -O ${results}/filtered_indels.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

# Part 3: Select only variants that have passed the filters
gatk SelectVariants \
    --exclude-filtered \
    -V ${results}/filtered_snps.vcf \
    -O ${results}/analysis-ready-snps.vcf

gatk SelectVariants \
    --exclude-filtered \
    -V ${results}/filtered_indels.vcf \
    -O ${results}/analysis-ready-indels.vcf

# Part 4: Remove variants that failed genotype filters for both SNPs and INDELs
grep -v -E "DP_filter|GQ_filter" ${results}/analysis-ready-snps.vcf > ${results}/analysis-ready-snps-filteredGT.vcf
grep -v -E "DP_filter|GQ_filter" ${results}/analysis-ready-indels.vcf > ${results}/analysis-ready-indels-filteredGT.vcf

# Part 5: Annotate filtered SNPs using GATK's Funcotator
gatk Funcotator \
    --variant ${results}/analysis-ready-snps-filteredGT.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path ${funcotator_data} \
    --output ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
    --output-file-format VCF

# Part 6: Annotate filtered INDELs using GATK's Funcotator
gatk Funcotator \
    --variant ${results}/analysis-ready-indels-filteredGT.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path ${funcotator_data} \
    --output ${results}/analysis-ready-indels-filteredGT-functotated.vcf \
    --output-file-format VCF

# Part 7: Export specific fields from SNPs annotated VCF to a tab-delimited table
gatk VariantsToTable \
    -V ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
    -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ${results}/output_snps.table
