# Germline Mutation Analysis Pipeline with GATK4 HaplotypeCaller

This repository contains a comprehensive pipeline script for **germline mutation analysis** using the Genome Analysis Toolkit (GATK4) HaplotypeCaller. The pipeline includes steps for variant calling, filtering, and annotation, demonstrating a complete bioinformatics workflow from downloading reference files and reads, to identifying and annotating germline mutations.

## Table of Contents
1. [Pipeline Overview](#pipeline-overview)
2. [Requirements](#requirements)
3. [Setup](#setup)
4. [Usage](#usage)
5. [Pipeline Steps](#pipeline-steps)
6. [Output](#output)
7. [Notes](#notes)

## Pipeline Overview
The pipeline is specifically designed for germline mutation calling in human or similar organisms. It includes:
- Downloading reference genome and known variants
- Aligning reads to a reference genome
- Marking duplicates and recalibrating base scores
- Calling germline variants (SNPs and INDELs) with GATK's HaplotypeCaller
- Filtering and annotating variants with Funcotator

This pipeline is ideal for researchers and bioinformaticians working on germline variant discovery.

## Requirements
This pipeline requires the following tools to be installed and accessible in your `$PATH`:
- **[bwa](http://bio-bwa.sourceforge.net/)** for sequence alignment
- **[samtools](http://www.htslib.org/)** for SAM/BAM file handling
- **[GATK4](https://gatk.broadinstitute.org/hc/en-us)** for germline variant calling, filtering, and annotation
- **wget** for downloading files (replaceable with other download tools)

To install these tools, consult their official documentation for installation instructions and dependencies.

## Setup
1. Clone this repository.
    ```bash
    git clone https://github.com/yourusername/germline-mutation-analysis-pipeline.git
    cd germline-mutation-analysis-pipeline
    ```

2. Ensure `bwa`, `samtools`, `GATK4`, and `wget` are installed.

3. Download any required data sources for GATK Funcotator annotation. Update paths in the script if needed.

## Usage
To run the pipeline:
1. Make the script executable:
    ```bash
    chmod +x pipeline.sh
    ```

2. Run the script:
    ```bash
    ./pipeline.sh
    ```

The script automatically downloads the necessary files (reference genome, known variants, and example reads) if they are not already present, so make sure you have an internet connection for the first run.

## Pipeline Steps
This script is divided into sections for clarity. Below is a breakdown of each part.

### PART 0: Download Required Files
The script downloads:
- Reference genome (FASTA file and associated index and dictionary files)
- Known sites for base quality score recalibration (BQSR) like dbSNP and known indels
- Sample input reads (paired-end FASTQ files for demonstration)
- GATK Funcotator data sources for annotation

### PART 1: Preprocessing and Germline Variant Calling
1. **Align Reads**: Uses `bwa mem` to align input reads to the reference genome.
2. **Convert and Sort**: Converts SAM to BAM, sorts, and indexes the BAM file.
3. **Mark Duplicates**: Identifies and marks duplicate reads using GATK.
4. **Base Quality Score Recalibration (BQSR)**: Recalibrates base scores based on known sites of variation.
5. **Germline Variant Calling**: Calls raw germline variants (SNPs and INDELs) using GATK HaplotypeCaller.

### PART 2: Filter and Annotate Variants
1. **Filter Variants**: Applies quality filters to SNPs and INDELs to reduce false positives.
2. **Select Passed Variants**: Extracts variants that passed filtering criteria for further analysis.
3. **Annotate with Funcotator**: Annotates filtered variants using GATK Funcotator.
4. **Export to Table**: Outputs key fields from annotated SNPs to a tab-delimited table.

## Output
Results are stored in the specified output directory (`output_directory/results`). Key output files include:
- `aligned_reads.sam` and `sorted_reads.bam`: Aligned read files
- `raw_variants.vcf`: Initial germline variant calls
- `filtered_snps.vcf` and `filtered_indels.vcf`: Quality-filtered SNPs and INDELs
- `analysis-ready-snps-filteredGT-functotated.vcf`: SNPs annotated with Funcotator
- `output_snps.table`: Tab-delimited table with annotated SNPs

## Notes
- **Adjust Paths**: Paths in the script (`ref`, `reads_dir`, `results`, `funcotator_data`) are customizable based on your file organization.
- **Download URLs**: Placeholder URLs for genome resources should be replaced with actual resource URLs (e.g., from [UCSC](http://genome.ucsc.edu/) or [NCBI](https://www.ncbi.nlm.nih.gov/)).
- **Runtime**: Processing times depend on data size and computational resources.
  
If you encounter any issues or have questions, feel free to open an issue or contact the repository maintainer.

---

This pipeline script is for educational and demonstrative purposes. Contributions to improve the script’s efficiency or add new features are welcome!
