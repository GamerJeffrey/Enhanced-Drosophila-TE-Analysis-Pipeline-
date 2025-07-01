# Enhanced-Drosophila-TE-Analysis-Pipeline-
# Enhanced Drosophila TE Analysis Pipeline

## Overview

This is a comprehensive, publication-ready pipeline for analyzing transposable element (TE) expression in single-cell RNA-seq data across multiple Drosophila species. The pipeline integrates cutting-edge bioinformatics tools with advanced statistical analysis to provide insights into TE expression patterns, cell-type specificity, and evolutionary conservation.

**✨ Remote Shell Optimized**: Designed to work on HPC clusters and remote systems without sudo access, with automatic detection of different GTF formats (NCBI vs FlyBase).

## Features

### Core Capabilities
- **Multi-species analysis**: Supports D. melanogaster, D. yakuba, and D. ananassae
- **Flexible file naming**: Automatically detects various GTF and FASTA naming conventions
- **Format awareness**: Handles both NCBI and FlyBase GTF formats automatically
- **No sudo required**: Works with conda, modules, or local installations
- **Comprehensive TE quantification**: Uses Telescope and featureCounts with consensus calling
- **Advanced clustering**: Multiple algorithms (Leiden, Louvain, hierarchical) with consensus
- **Quality control**: Extensive QC with FastQC, cell filtering, and doublet detection
- **Cross-species comparison**: Orthology-based comparative analysis
- **Evolutionary insights**: Phylogenetic context and species-specific patterns

### Remote System Features
- **HPC compatibility**: Works with module systems and job schedulers
- **Fallback mechanisms**: Graceful degradation when tools are missing
- **Local installation**: Can install tools in user directory
- **Memory efficient**: Configurable resource usage
- **Progress monitoring**: Comprehensive logging and progress tracking

## Quick Start for Remote Systems

### 1. Setup Environment (First Time Only)
```bash
# Download and setup the pipeline
git clone [repository-url] te_pipeline
cd te_pipeline

# Setup tools locally (if needed)
bash setup_remote_env.sh

# Install R packages
Rscript install_r_packages.R

# Activate environment
source $HOME/te_pipeline/activate.sh
```

### 2. Prepare Your Data
```bash
# Create directory structure
mkdir -p references fastqs

# Your file structure can be flexible:
# references/
#   dmel_genome.fasta, dmel_genes.gtf, dmel_TEs.gtf
#   dyak_genome.fasta, dyak_genes.gtf, dyak_repeats.gtf  
#   dana_genome.fasta, dana_genes.gtf, dana_repeats.gtf
# fastqs/
#   Dmel/sample_R1.fastq.gz, sample_R2.fastq.gz
#   Dyak/sample_R1.fastq.gz, sample_R2.fastq.gz
#   Dana/sample_R1.fastq.gz, sample_R2.fastq.gz
```

### 3. Run Pipeline
```bash
# Basic run
./enhanced_te_pipeline.sh

# With custom config
./enhanced_te_pipeline.sh config/my_config.yaml

# Check what tools are available
./enhanced_te_pipeline.sh --check-deps
```

### Advanced Analysis
- **Trajectory analysis**: Pseudotime analysis of TE expression dynamics
- **Co-expression networks**: TE-gene co-expression analysis
- **Differential expression**: Statistical testing with multiple correction methods
- **Publication-ready figures**: Automated generation of high-quality visualizations
- **Comprehensive reporting**: HTML reports with interactive elements

### Quality & Reproducibility
- **Resource monitoring**: CPU and memory usage tracking
- **Version control**: Complete software version logging
- **Error handling**: Robust error checking and cleanup
- **Parallel processing**: Efficient multi-core utilization
- **Validation**: Comprehensive result validation
- **Checkpointing**: Resume from intermediate steps

## Installation for Remote Systems

### Prerequisites
- Linux or macOS system
- Access to compute resources (8+ GB RAM recommended)
- Internet access for downloading tools and databases

### Installation Options

#### Option 1: Conda (Recommended)
```bash
# Install Miniconda if not available
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

# Create and activate environment
conda create -n te_pipeline python=3.9
conda activate te_pipeline

# Install tools
conda install -c bioconda -c conda-forge star samtools subread parallel
pip install telescope-bth
```

#### Option 2: HPC Module System
```bash
# Load modules (adjust names for your system)
module load star/2.7.10 samtools/1.15 subread/2.0.3 R/4.2.0 python/3.9

# Add to your job scripts
echo "module load star samtools subread R python" >> load_modules.sh
```

#### Option 3: Local Installation
```bash
# Use the provided setup script
bash setup_remote_env.sh
```

### Supported File Naming Patterns

The pipeline automatically detects various naming conventions:

#### Genome Files
- `{species}_genome.fasta` (e.g., `dmel_genome.fasta`)

#### Gene Annotation Files  
- `{species}_genes.gtf`
- `{species}.gtf` 
- `genes_{species}.gtf`
- Any file matching `*{species}*genes*.gtf`

#### TE Annotation Files
- `{species}TEs.gtf`
- `{species}_TEs.gtf`
- `{species}_repeats.gtf`
- Any file matching `*{species}*TE*.gtf` or `*{species}*repeat*.gtf`

#### FASTQ Files
- Located in `fastqs/{Species}/` directories
- Paired-end: `*_R1.fastq.gz` and `*_R2.fastq.gz`

### GTF Format Compatibility

The pipeline automatically handles different GTF formats:

- **FlyBase format** (D. melanogaster): Standard FlyBase attributes
- **NCBI format** (D. yakuba, D. ananassae): NCBI RefSeq attributes
- **Custom formats**: Automatically converts attributes as needed
