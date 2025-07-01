#SBATCH --job-name=Enhanced_TE_Pipe
#SBATCH --output=Enhanced_TE_Pipe_%j.out
#SBATCH --error=Enhanced_TE_Pipe_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general

# Enhanced Drosophila TE Analysis Pipeline
# Publication-ready version with comprehensive analysis and quality control
# Author: Enhanced Pipeline v2.0
# Date: $(date)

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Global variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_VERSION="2.0"
START_TIME=$(date +%s)

# Default configuration
DEFAULT_CONFIG="config/pipeline_config.yaml"
LOG_DIR="logs"
TEMP_DIR="temp"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_DIR/pipeline.log"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_DIR/pipeline.log"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_DIR/pipeline.log"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_DIR/pipeline.log"
}

# Enhanced cleanup function
cleanup_on_error() {
    local exit_code=$?
    log_error "Pipeline failed with exit code $exit_code. Performing cleanup..."
    
    # Kill background processes
    jobs -p | xargs -r kill 2>/dev/null || true
    
    # Clean temporary files
    if [ -d "$TEMP_DIR" ]; then
        rm -rf "$TEMP_DIR"
    fi
    
    # Remove incomplete output files
    find "$OUTDIR" -name "*.tmp" -delete 2>/dev/null || true
    find "$OUTDIR" -name "*.partial" -delete 2>/dev/null || true
    
    log_error "Cleanup completed. Check logs for details."
    exit $exit_code
}

trap cleanup_on_error ERR INT TERM

# Configuration loading
load_config() {
    local config_file=${1:-$DEFAULT_CONFIG}
    
    if [ ! -f "$config_file" ]; then
        log_warning "Config file $config_file not found. Creating default configuration."
        create_default_config "$config_file"
    fi
    
    log_info "Loading configuration from $config_file"
    
    # Parse YAML-like config (simplified)
    eval $(grep -E '^[[:space:]]*[^#].*:' "$config_file" | \
           sed 's/[[:space:]]*//g' | \
           sed 's/:/=/g' | \
           sed 's/^/export /')
    
    # Set defaults if not specified
    export THREADS=${THREADS:-8}
    export MEMORY=${MEMORY:-32G}
    export MIN_CELLS=${MIN_CELLS:-10}
    export MIN_GENES=${MIN_GENES:-200}
    export MAX_MITO_PERCENT=${MAX_MITO_PERCENT:-20}
    export RESOLUTION=${RESOLUTION:-0.5}
    export MIN_TE_COUNT=${MIN_TE_COUNT:-5}
    export SEED=${SEED:-42}
    export SPECIES_LIST=${SPECIES_LIST:-"Dmel,Dyak,Dana"}
    
    # Convert comma-separated to array
    IFS=',' read -ra SPECIES_ARRAY <<< "$SPECIES_LIST"
    
    log_success "Configuration loaded successfully"
}

# Create default configuration file
create_default_config() {
    local config_file=$1
    mkdir -p "$(dirname "$config_file")"
    
    cat > "$config_file" << 'EOF'
# Enhanced Drosophila TE Analysis Pipeline Configuration

# Computational resources
threads: 8
memory: 32G

# Quality control parameters
min_cells: 10
min_genes: 200
max_mito_percent: 20

# Clustering parameters
resolution: 0.5
seed: 42

# TE analysis parameters
min_te_count: 5
te_families_of_interest: "LINE,LTR,DNA,SINE"

# Species to analyze
species_list: "Dmel,Dyak,Dana"

# Analysis options
run_trajectory: true
run_coexpression: true
run_differential: true
run_evolutionary: true

# Output options
create_html_report: true
compress_outputs: true

# External database URLs
flybase_orthology_url: "ftp://ftp.flybase.net/releases/FB2023_05/precomputed_files/orthologs/dmel_orthologs_in_drosophila_species_fb_2023_05.tsv.gz"
dfam_url: "https://dfam.org/releases/Dfam_3.7/families/Dfam_3.7_families.tsv.gz"
EOF
    
    log_success "Default configuration created at $config_file"
}

# System requirements check without sudo/brew dependencies
check_system_requirements() {
    log_info "Checking system requirements..."
    
    # Check for required tools that can be installed without sudo
    local required_tools=("samtools" "Rscript" "python3" "wget" "curl")
    local optional_tools=("STAR" "telescope" "featureCounts" "parallel")
    local missing_required=()
    local missing_optional=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_required+=("$tool")
        fi
    done
    
    for tool in "${optional_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_optional+=("$tool")
        fi
    done
    
    if [ ${#missing_required[@]} -ne 0 ]; then
        log_error "Missing required tools: ${missing_required[*]}"
        log_error "Please install missing tools using conda or module load, then try again."
        log_info "Example: conda install -c bioconda ${missing_required[*]}"
        exit 1
    fi
    
    if [ ${#missing_optional[@]} -ne 0 ]; then
        log_warning "Missing optional tools: ${missing_optional[*]}"
        log_warning "Some functionality may be limited. Consider installing with conda:"
        log_warning "conda install -c bioconda ${missing_optional[*]}"
    fi
    
    # Check STAR specifically and suggest alternatives
    if ! command -v STAR &> /dev/null; then
        log_warning "STAR not found. Checking for alternatives..."
        if command -v module &> /dev/null; then
            log_info "Try: module load star  # or similar module name"
        fi
        log_info "Or install with: conda install -c bioconda star"
    fi
    
    # Check available memory
    local available_memory=""
    if command -v free &> /dev/null; then
        available_memory=$(free -g | awk '/^Mem:/{print $7}')
    elif [ -f /proc/meminfo ]; then
        available_memory=$(awk '/MemAvailable/ {print int($2/1024/1024)}' /proc/meminfo)
    fi
    
    if [ -n "$available_memory" ]; then
        local required_memory=${MEMORY%G}
        if [ "$available_memory" -lt "$required_memory" ]; then
            log_warning "Available memory ($available_memory GB) is less than requested ($required_memory GB)"
            log_warning "Consider reducing memory usage in config: memory: ${available_memory}G"
        fi
    fi
    
    # Check disk space
    local available_space=""
    if command -v df &> /dev/null; then
        available_space=$(df -BG "$PWD" 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//' || echo "unknown")
    fi
    
    if [ "$available_space" != "unknown" ] && [ "$available_space" -lt 100 ]; then
        log_warning "Available disk space is less than 100GB. Pipeline may fail due to insufficient storage."
    fi
    
    # Check if we can create temporary files
    if ! touch "$TEMP_DIR/test_write" 2>/dev/null; then
        log_error "Cannot write to temporary directory: $TEMP_DIR"
        exit 1
    fi
    rm -f "$TEMP_DIR/test_write"
    
    log_success "System requirements check completed"
}

# Enhanced file checking with detailed validation and flexible naming
check_files() {
    local species=$1
    log_info "Validating input files for $species..."
    
    # Flexible file naming patterns
    local genome_file="$REFDIR/${species,,}_genome.fasta"
    local genes_file=""
    local te_file=""
    local fastq_dir="$FASTQDIR/$species"
    
    # Find genes file with flexible naming
    if [ -f "$REFDIR/${species,,}_genes.harmonized.gtf" ]; then
        genes_file="$REFDIR/${species,,}_genes.harmonized.gtf"
    elif [ -f "$REFDIR/${species,,}_genes.gtf" ]; then
        genes_file="$REFDIR/${species,,}_genes.gtf"
    elif [ -f "$REFDIR/${species,,}.gtf" ]; then
        genes_file="$REFDIR/${species,,}.gtf"
    elif [ -f "$REFDIR/genes_${species,,}.gtf" ]; then
        genes_file="$REFDIR/genes_${species,,}.gtf"
    else
        # Search for any GTF file containing the species name
        genes_file=$(find "$REFDIR" -name "*${species,,}*genes*.gtf" -o -name "*${species,,}*.gtf" | head -1)
    fi
    
    # Find TE file with flexible naming
    if [ -f "$REFDIR/${species,,}TEs.harmonized.gtf" ]; then
        te_file="$REFDIR/${species,,}TEs.harmonized.gtf"
    elif [ -f "$REFDIR/${species,,}_TEs.gtf" ]; then
        te_file="$REFDIR/${species,,}_TEs.gtf"
    elif [ -f "$REFDIR/${species,,}_repeats.gtf" ]; then
        te_file="$REFDIR/${species,,}_repeats.gtf"
    elif [ -f "$REFDIR/TEs_${species,,}.gtf" ]; then
        te_file="$REFDIR/TEs_${species,,}.gtf"
    else
        # Search for any GTF file containing TE/repeat keywords
        te_file=$(find "$REFDIR" -name "*${species,,}*TE*.gtf" -o -name "*${species,,}*repeat*.gtf" | head -1)
    fi
    
    log_info "Found files for $species:"
    log_info "  Genome: $genome_file"
    log_info "  Genes: $genes_file"
    log_info "  TEs: $te_file"
    
    local errors=0
    
    # Check genome file
    if [ ! -f "$genome_file" ]; then
        log_error "Missing genome file: $genome_file"
        ((errors++))
    else
        # Validate FASTA format
        if ! head -1 "$genome_file" | grep -q "^>"; then
            log_error "Invalid FASTA format in $genome_file"
            ((errors++))
        else
            local seq_count=$(grep -c "^>" "$genome_file")
            log_info "Genome file contains $seq_count sequences"
        fi
    fi
    
    # Check GTF files with format detection
    for gtf_file in "$genes_file" "$te_file"; do
        if [ -z "$gtf_file" ] || [ ! -f "$gtf_file" ]; then
            log_error "Missing GTF file: $gtf_file"
            ((errors++))
        else
            # Detect GTF format (NCBI vs FlyBase)
            local gtf_format="unknown"
            if grep -q "gene_id" "$gtf_file" | head -1; then
                if grep -q "FBgn" "$gtf_file" | head -1; then
                    gtf_format="flybase"
                else
                    gtf_format="ncbi"
                fi
            fi
            
            # Validate GTF format
            if ! head -1 "$gtf_file" | grep -qE 
    
    # Check FASTQ directory and files
    if [ ! -d "$fastq_dir" ]; then
        log_error "Missing FASTQ directory: $fastq_dir"
        ((errors++))
    else
        local fastq_files=($(find "$fastq_dir" -name "*.fastq.gz" | sort))
        if [ ${#fastq_files[@]} -lt 2 ]; then
            log_error "Expected at least 2 FASTQ files in $fastq_dir, found ${#fastq_files[@]}"
            ((errors++))
        else
            log_info "Found ${#fastq_files[@]} FASTQ files for $species"
            
            # Validate FASTQ format and check for paired-end
            for fq in "${fastq_files[@]:0:2}"; do
                if ! zcat "$fq" 2>/dev/null | head -4 | grep -q "^@"; then
                    log_error "Invalid FASTQ format in $fq"
                    ((errors++))
                fi
            done
        fi
    fi
    
    if [ $errors -eq 0 ]; then
        log_success "All files validated for $species"
        return 0
    else
        log_error "$errors validation errors found for $species"
        return 1
    fi
}

# Resource monitoring function
monitor_resources() {
    local pid=$1
    local process_name=$2
    local logfile="$LOG_DIR/resources_${process_name}.log"
    
    echo "Time,PID,CPU%,MEM%,VSZ,RSS" > "$logfile"
    
    while kill -0 $pid 2>/dev/null; do
        ps -p $pid -o pid,pcpu,pmem,vsz,rss --no-headers >> "$logfile" 2>/dev/null || break
        sleep 30
    done &
    
    echo $!  # Return monitoring PID
}

# Enhanced quality control with FastQC
run_fastqc() {
    local species=$1
    log_info "Running FastQC quality control for $species..."
    
    local fastq_dir="$FASTQDIR/$species"
    local qc_dir="$OUTDIR/qc/fastqc_${species}"
    
    mkdir -p "$qc_dir"
    
    # Run FastQC on all FASTQ files
    fastqc -t "$THREADS" -o "$qc_dir" "$fastq_dir"/*.fastq.gz
    
    # Generate summary report
    python3 "$SCRIPTSDIR/summarize_fastqc.py" \
        --input_dir "$qc_dir" \
        --output "$qc_dir/fastqc_summary.json"
    
    log_success "FastQC completed for $species"
}

# Version tracking
record_versions() {
    local version_file="$OUTDIR/software_versions.txt"
    
    log_info "Recording software versions..."
    
    cat > "$version_file" << EOF
# Software Versions Used in Analysis
# Generated on: $(date)
# Pipeline Version: $PIPELINE_VERSION

EOF
    
    # Record versions of key tools
    echo "STAR: $(STAR --version 2>&1 | head -1)" >> "$version_file"
    echo "samtools: $(samtools --version | head -1)" >> "$version_file"
    echo "R: $(Rscript --version 2>&1 | head -1)" >> "$version_file"
    echo "Python: $(python3 --version)" >> "$version_file"
    
    # Record R package versions
    Rscript -e "
    packages <- c('Seurat', 'dplyr', 'ggplot2', 'SingleCellExperiment', 'scater', 'SAKE', 'DESeq2', 'monocle3')
    for(pkg in packages) {
        if(require(pkg, quietly=TRUE, character.only=TRUE)) {
            cat(paste0(pkg, ': ', packageVersion(pkg), '\n'))
        }
    }
    " >> "$version_file"
    
    log_success "Software versions recorded"
}

# Enhanced STAR index building with format-aware GTF handling
build_star_index() {
    local species=$1
    log_info "Building STAR index for $species..."
    
    local genome_file="$REFDIR/${species,,}_genome.fasta"
    local genes_file=$(find "$REFDIR" -name "*${species,,}*genes*.gtf" -o -name "*${species,,}*.gtf" | grep -v TE | grep -v repeat | head -1)
    local index_dir="$OUTDIR/star_indices/star_index_${species}"
    
    if [ -z "$genes_file" ]; then
        log_error "Could not find genes GTF file for $species"
        return 1
    fi
    
    mkdir -p "$index_dir"
    
    # Detect GTF format and apply appropriate processing
    local gtf_format=$(grep "$species:genes:" "$TEMP_DIR/gtf_formats.txt" 2>/dev/null | cut -d: -f3)
    
    # Create processed GTF if needed
    local processed_gtf="$TEMP_DIR/${species}_genes_processed.gtf"
    
    if [ "$gtf_format" = "ncbi" ]; then
        log_info "Processing NCBI format GTF for $species..."
        # Convert NCBI GTF to be compatible with STAR
        awk -F'\t' 'BEGIN{OFS="\t"} 
        $3=="exon" {
            # Extract gene_id and transcript_id from attributes
            if (match($9, /gene_id "([^"]+)"/, gene_arr)) {
                gene_id = gene_arr[1]
            }
            if (match($9, /transcript_id "([^"]+)"/, trans_arr)) {
                transcript_id = trans_arr[1]
            } else {
                transcript_id = gene_id "_1"
            }
            # Rebuild attributes in standard format
            $9 = "gene_id \"" gene_id "\"; transcript_id \"" transcript_id "\";"
            print $0
        }' "$genes_file" > "$processed_gtf"
    else
        # Use original GTF for FlyBase format
        cp "$genes_file" "$processed_gtf"
    fi
    
    # Calculate optimal parameters based on genome size
    local genome_size=$(stat -c%s "$genome_file" 2>/dev/null || stat -f%z "$genome_file" 2>/dev/null || echo "2000000000")
    local genomeSAindexNbases=$((14))  # Default for most genomes
    
    if [ $genome_size -lt 2000000000 ]; then  # < 2GB
        genomeSAindexNbases=13
    fi
    
    # Run STAR with resource monitoring
    STAR --runMode genomeGenerate \
         --genomeDir "$index_dir" \
         --genomeFastaFiles "$genome_file" \
         --sjdbGTFfile "$processed_gtf" \
         --sjdbOverhang 100 \
         --runThreadN "$THREADS" \
         --limitGenomeGenerateRAM $(($(echo $MEMORY | sed 's/G//') * 1000000000)) \
         --genomeSAindexNbases $genomeSAindexNbases \
         --outTmpDir "$TEMP_DIR/star_tmp_${species}" &
    
    local star_pid=$!
    local monitor_pid=$(monitor_resources $star_pid "star_index_${species}")
    
    wait $star_pid
    kill $monitor_pid 2>/dev/null || true
    
    # Validate index
    if [ ! -f "$index_dir/SA" ]; then
        log_error "STAR index building failed for $species"
        return 1
    fi
    
    log_success "STAR index built successfully for $species"
}

# Enhanced STARsolo with comprehensive parameters
run_starsolo() {
    local species=$1
    log_info "Running STARsolo alignment for $species..."
    
    local fastq_dir="$FASTQDIR/$species"
    local fastq_files=($(find "$fastq_dir" -name "*.fastq.gz" | sort))
    local fq1="${fastq_files[0]}"
    local fq2="${fastq_files[1]}"
    local index_dir="$OUTDIR/star_indices/star_index_${species}"
    local output_dir="$OUTDIR/alignments/starsolo_${species}"
    
    mkdir -p "$output_dir"
    
    # Enhanced STARsolo parameters
    STAR --genomeDir "$index_dir" \
         --readFilesIn "$fq1" "$fq2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$output_dir/" \
         --soloType CB_UMI_Simple \
         --soloCBwhitelist None \
         --soloCBstart 1 --soloCBlen 16 \
         --soloUMIstart 17 --soloUMIlen 12 \
         --soloStrand Forward \
         --soloFeatures Gene GeneFull Velocyto \
         --soloOutFileNames Solo.out/ features.tsv barcodes.tsv matrix.mtx \
         --soloMultiMappers EM \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --runThreadN "$THREADS" \
         --limitBAMsortRAM $(($(echo $MEMORY | sed 's/G//') * 1000000000)) \
         --outTmpDir "$TEMP_DIR/starsolo_tmp_${species}" &
    
    local star_pid=$!
    local monitor_pid=$(monitor_resources $star_pid "starsolo_${species}")
    
    wait $star_pid
    kill $monitor_pid 2>/dev/null || true
    
    # Copy and index the sorted BAM file
    cp "$output_dir/Aligned.sortedByCoord.out.bam" "$OUTDIR/alignments/starsolo_${species}.bam"
    samtools index "$OUTDIR/alignments/starsolo_${species}.bam"
    
    # Generate alignment statistics
    collect_alignment_stats "$species"
    
    log_success "STARsolo alignment completed for $species"
}

# Comprehensive alignment statistics
collect_alignment_stats() {
    local species=$1
    log_info "Collecting alignment statistics for $species..."
    
    local bam_file="$OUTDIR/alignments/starsolo_${species}.bam"
    local stats_dir="$OUTDIR/qc/alignment_stats"
    mkdir -p "$stats_dir"
    
    # Basic alignment statistics
    samtools flagstat "$bam_file" > "$stats_dir/${species}_flagstat.txt"
    samtools idxstats "$bam_file" > "$stats_dir/${species}_idxstats.txt"
    samtools stats "$bam_file" > "$stats_dir/${species}_stats.txt"
    
    # STARsolo statistics
    local starsolo_dir="$OUTDIR/alignments/starsolo_${species}"
    if [ -f "$starsolo_dir/Log.final.out" ]; then
        cp "$starsolo_dir/Log.final.out" "$stats_dir/${species}_star_final.log"
    fi
    
    # Generate summary statistics
    python3 "$SCRIPTSDIR/summarize_alignment_stats.py" \
        --species "$species" \
        --stats_dir "$stats_dir" \
        --output "$stats_dir/${species}_alignment_summary.json"
    
    log_success "Alignment statistics collected for $species"
}

# Enhanced TE quantification with format-aware processing
run_te_quantification() {
    local species=$1
    log_info "Running comprehensive TE quantification for $species..."
    
    local bam_file="$OUTDIR/alignments/starsolo_${species}.bam"
    local te_file=$(find "$REFDIR" -name "*${species,,}*TE*.gtf" -o -name "*${species,,}*repeat*.gtf" | head -1)
    local output_dir="$OUTDIR/te_counts"
    
    if [ -z "$te_file" ]; then
        log_error "Could not find TE GTF file for $species"
        return 1
    fi
    
    mkdir -p "$output_dir"
    
    # Process TE GTF based on format
    local gtf_format=$(grep "$species:tes:" "$TEMP_DIR/gtf_formats.txt" 2>/dev/null | cut -d: -f3)
    local processed_te_gtf="$TEMP_DIR/${species}_tes_processed.gtf"
    
    if [ "$gtf_format" = "ncbi" ]; then
        log_info "Processing NCBI format TE GTF for $species..."
        # Convert NCBI TE GTF format
        awk -F'\t' 'BEGIN{OFS="\t"} 
        $3=="exon" || $3=="repeat_region" {
            # Extract or create gene_id and transcript_id
            if (match($9, /gene_id "([^"]+)"/, gene_arr)) {
                gene_id = gene_arr[1]
            } else if (match($9, /Name=([^;]+)/, name_arr)) {
                gene_id = name_arr[1]
            } else {
                gene_id = $1 "_" $4 "_" $5
            }
            
            if (match($9, /transcript_id "([^"]+)"/, trans_arr)) {
                transcript_id = trans_arr[1]
            } else {
                transcript_id = gene_id "_1"
            }
            
            # Convert to exon feature for counting
            $3 = "exon"
            $9 = "gene_id \"" gene_id "\"; transcript_id \"" transcript_id "\";"
            print $0
        }' "$te_file" > "$processed_te_gtf"
    else
        # Use original GTF for FlyBase format
        cp "$te_file" "$processed_te_gtf"
    fi
    
    # Sort BAM by name for Telescope
    local name_sorted_bam="$OUTDIR/alignments/starsolo_${species}.sortedByName.bam"
    samtools sort -n -@ "$THREADS" -o "$name_sorted_bam" "$bam_file"
    
    # Run Telescope (if available)
    if command -v telescope &> /dev/null; then
        log_info "Running Telescope for $species..."
        telescope assign \
            --barcode_tag CB \
            --umi_tag UB \
            --stranded_library \
            --theta_prior 200000 \
            --max_iter 200 \
            --updated_sam \
            "$name_sorted_bam" \
            "$processed_te_gtf" \
            > "$output_dir/telescope_${species}_report.tsv" 2>/dev/null &
        
        local telescope_pid=$!
    else
        log_warning "Telescope not found, skipping Telescope quantification"
        local telescope_pid=""
    fi
    
    # Run featureCounts for comparison
    if command -v featureCounts &> /dev/null; then
        log_info "Running featureCounts for $species..."
        featureCounts \
            -T "$THREADS" \
            -p \
            -B \
            -C \
            -s 1 \
            -a "$processed_te_gtf" \
            -o "$output_dir/featurecounts_${species}_te.txt" \
            "$bam_file" 2>/dev/null &
        
        local fc_pid=$!
    else
        log_warning "featureCounts not found, skipping featureCounts quantification"
        local fc_pid=""
    fi
    
    # Wait for both methods
    if [ -n "$telescope_pid" ]; then
        wait $telescope_pid
    fi
    if [ -n "$fc_pid" ]; then
        wait $fc_pid
    fi
    
    # Process and harmonize results
    Rscript "$SCRIPTSDIR/harmonize_te_counts.R" \
        --telescope "$output_dir/telescope_${species}_report.tsv" \
        --featurecounts "$output_dir/featurecounts_${species}_te.txt" \
        --species "$species" \
        --output "$output_dir/harmonized_te_counts_${species}.rds"
    
    log_success "TE quantification completed for $species"
}

# Advanced cell filtering and quality control
run_cell_qc() {
    local species=$1
    log_info "Running advanced cell quality control for $species..."
    
    local starsolo_dir="$OUTDIR/alignments/starsolo_${species}"
    local qc_dir="$OUTDIR/qc/cell_qc_${species}"
    
    mkdir -p "$qc_dir"
    
    Rscript "$SCRIPTSDIR/advanced_cell_qc.R" \
        --input "$starsolo_dir/Solo.out/Gene/filtered/" \
        --species "$species" \
        --min_cells "$MIN_CELLS" \
        --min_genes "$MIN_GENES" \
        --max_mito_percent "$MAX_MITO_PERCENT" \
        --output "$qc_dir/" \
        --seed "$SEED"
    
    log_success "Cell quality control completed for $species"
}

# Enhanced clustering with multiple algorithms
run_advanced_clustering() {
    local species=$1
    log_info "Running advanced clustering analysis for $species..."
    
    local qc_dir="$OUTDIR/qc/cell_qc_${species}"
    local cluster_dir="$OUTDIR/clusters/advanced_${species}"
    
    mkdir -p "$cluster_dir"
    
    # Run multiple clustering algorithms
    Rscript "$SCRIPTSDIR/advanced_clustering.R" \
        --input "$qc_dir/filtered_data.rds" \
        --species "$species" \
        --resolution "$RESOLUTION" \
        --algorithms "leiden,louvain,hierarchical" \
        --output "$cluster_dir/" \
        --seed "$SEED"
    
    log_success "Advanced clustering completed for $species"
}

# Trajectory analysis
run_trajectory_analysis() {
    local species=$1
    
    if [ "$run_trajectory" != "true" ]; then
        log_info "Skipping trajectory analysis for $species (disabled in config)"
        return 0
    fi
    
    log_info "Running trajectory analysis for $species..."
    
    local cluster_dir="$OUTDIR/clusters/advanced_${species}"
    local trajectory_dir="$OUTDIR/trajectory/${species}"
    
    mkdir -p "$trajectory_dir"
    
    Rscript "$SCRIPTSDIR/trajectory_analysis.R" \
        --input "$cluster_dir/clustered_data.rds" \
        --species "$species" \
        --output "$trajectory_dir/" \
        --seed "$SEED"
    
    log_success "Trajectory analysis completed for $species"
}

# Co-expression network analysis
run_coexpression_analysis() {
    local species=$1
    
    if [ "$run_coexpression" != "true" ]; then
        log_info "Skipping co-expression analysis for $species (disabled in config)"
        return 0
    fi
    
    log_info "Running co-expression network analysis for $species..."
    
    local te_counts="$OUTDIR/te_counts/harmonized_te_counts_${species}.rds"
    local gene_counts="$OUTDIR/alignments/starsolo_${species}/Solo.out/Gene/filtered/"
    local coexp_dir="$OUTDIR/coexpression/${species}"
    
    mkdir -p "$coexp_dir"
    
    Rscript "$SCRIPTSDIR/coexpression_analysis.R" \
        --te_counts "$te_counts" \
        --gene_counts "$gene_counts" \
        --species "$species" \
        --output "$coexp_dir/" \
        --seed "$SEED"
    
    log_success "Co-expression analysis completed for $species"
}

# Enhanced orthology analysis with phylogenetic context
setup_enhanced_orthology() {
    log_info "Setting up enhanced orthology analysis..."
    
    local ortho_dir="$REFDIR/orthology"
    mkdir -p "$ortho_dir"
    cd "$ortho_dir"
    
    # Download FlyBase orthology data
    if [ ! -f "dmel_orthologs_in_drosophila_species_fb_2023_05.tsv" ]; then
        log_info "Downloading FlyBase orthology data..."
        wget -O dmel_orthologs_in_drosophila_species_fb_2023_05.tsv.gz \
            "$flybase_orthology_url"
        gunzip dmel_orthologs_in_drosophila_species_fb_2023_05.tsv.gz
    fi
    
    # Download Dfam data for TE classification
    if [ ! -f "Dfam_3.7_families.tsv" ]; then
        log_info "Downloading Dfam TE classification data..."
        wget -O Dfam_3.7_families.tsv.gz "$dfam_url"
        gunzip Dfam_3.7_families.tsv.gz
    fi
    
    # Download phylogenetic tree
    if [ ! -f "drosophila_species_tree.nwk" ]; then
        log_info "Creating phylogenetic tree for Drosophila species..."
        cat > drosophila_species_tree.nwk << 'EOF'
(Dmel:0.1,(Dyak:0.05,Dana:0.05):0.05);
EOF
    fi
    
    cd - > /dev/null
    log_success "Enhanced orthology setup completed"
}

# Comprehensive orthology processing with evolutionary analysis
process_enhanced_orthology() {
    log_info "Processing enhanced orthology data with evolutionary analysis..."
    
    Rscript "$SCRIPTSDIR/enhanced_orthology_analysis.R" \
        --orthology_file "$REFDIR/orthology/dmel_orthologs_in_drosophila_species_fb_2023_05.tsv" \
        --dfam_file "$REFDIR/orthology/Dfam_3.7_families.tsv" \
        --tree_file "$REFDIR/orthology/drosophila_species_tree.nwk" \
        --species_list "${SPECIES_ARRAY[*]}" \
        --output "$OUTDIR/orthology/"
    
    log_success "Enhanced orthology processing completed"
}

# Comprehensive differential analysis
run_differential_analysis() {
    if [ "$run_differential" != "true" ]; then
        log_info "Skipping differential analysis (disabled in config)"
        return 0
    fi
    
    log_info "Running comprehensive differential analysis..."
    
    local diff_dir="$OUTDIR/differential_analysis"
    mkdir -p "$diff_dir"
    
    # Collect all integration files
    local integration_files=()
    for species in "${SPECIES_ARRAY[@]}"; do
        integration_files+=("$OUTDIR/integration_${species}.rds")
    done
    
    Rscript "$SCRIPTSDIR/comprehensive_differential_analysis.R" \
        --integration_files "${integration_files[@]}" \
        --orthology_data "$OUTDIR/orthology/processed_orthology.rds" \
        --output "$diff_dir/" \
        --species_list "${SPECIES_ARRAY[*]}" \
        --seed "$SEED"
    
    log_success "Differential analysis completed"
}

# Evolutionary analysis
run_evolutionary_analysis() {
    if [ "$run_evolutionary" != "true" ]; then
        log_info "Skipping evolutionary analysis (disabled in config)"
        return 0
    fi
    
    log_info "Running evolutionary analysis..."
    
    local evo_dir="$OUTDIR/evolutionary_analysis"
    mkdir -p "$evo_dir"
    
    Rscript "$SCRIPTSDIR/evolutionary_analysis.R" \
        --integration_dir "$OUTDIR/" \
        --orthology_data "$OUTDIR/orthology/processed_orthology.rds" \
        --tree_file "$REFDIR/orthology/drosophila_species_tree.nwk" \
        --output "$evo_dir/" \
        --species_list "${SPECIES_ARRAY[*]}"
    
    log_success "Evolutionary analysis completed"
}

# Enhanced integration of TE counts with clusters
integrate_te_clusters_enhanced() {
    local species=$1
    log_info "Running enhanced TE-cluster integration for $species..."
    
    local te_counts="$OUTDIR/te_counts/harmonized_te_counts_${species}.rds"
    local cluster_data="$OUTDIR/clusters/advanced_${species}/clustered_data.rds"
    local cell_qc="$OUTDIR/qc/cell_qc_${species}/cell_metadata.rds"
    local output_file="$OUTDIR/integration_${species}.rds"
    
    Rscript "$SCRIPTSDIR/enhanced_te_cluster_integration.R" \
        --te_counts "$te_counts" \
        --cluster_data "$cluster_data" \
        --cell_qc "$cell_qc" \
        --species "$species" \
        --output "$output_file" \
        --min_te_count "$MIN_TE_COUNT"
    
    log_success "Enhanced TE-cluster integration completed for $species"
}

# Publication-ready figure generation
create_publication_figures() {
    log_info "Creating publication-ready figures..."
    
    local figures_dir="$OUTDIR/figures"
    mkdir -p "$figures_dir"
    
    # Collect all integration files
    local integration_files=()
    for species in "${SPECIES_ARRAY[@]}"; do
        if [ -f "$OUTDIR/integration_${species}.rds" ]; then
            integration_files+=("$OUTDIR/integration_${species}.rds")
        fi
    done
    
    Rscript "$SCRIPTSDIR/create_publication_figures.R" \
        --integration_files "${integration_files[@]}" \
        --orthology_data "$OUTDIR/orthology/processed_orthology.rds" \
        --differential_data "$OUTDIR/differential_analysis/" \
        --evolutionary_data "$OUTDIR/evolutionary_analysis/" \
        --output "$figures_dir/" \
        --species_list "${SPECIES_ARRAY[*]}"
    
    log_success "Publication-ready figures created"
}

# Comprehensive HTML report generation
generate_html_report() {
    if [ "$create_html_report" != "true" ]; then
        log_info "Skipping HTML report generation (disabled in config)"
        return 0
    fi
    
    log_info "Generating comprehensive HTML report..."
    
    local report_dir="$OUTDIR/reports"
    mkdir -p "$report_dir"
    
    Rscript "$SCRIPTSDIR/generate_html_report.R" \
        --pipeline_dir "$OUTDIR" \
        --config_file "$DEFAULT_CONFIG" \
        --versions_file "$OUTDIR/software_versions.txt" \
        --output "$report_dir/pipeline_report.html" \
        --species_list "${SPECIES_ARRAY[*]}"
    
    log_success "HTML report generated at $report_dir/pipeline_report.html"
}

# Fallback parallel processing for systems without GNU parallel
run_parallel_analysis_fallback() {
    local func_name=$1
    shift
    local species_list=("$@")
    
    log_info "Running $func_name sequentially for ${#species_list[@]} species (parallel not available)"
    
    for species in "${species_list[@]}"; do
        log_info "Processing $species with $func_name..."
        eval "$func_name $species"
    done
    
    log_success "Sequential $func_name completed for all species"
}

# Enhanced parallel processing wrapper with fallback
run_parallel_analysis() {
    local func_name=$1
    shift
    local species_list=("$@")
    
    if command -v parallel &> /dev/null; then
        log_info "Running $func_name in parallel for ${#species_list[@]} species"
        
        # Create array of commands
        local commands=()
        for species in "${species_list[@]}"; do
            commands+=("$func_name $species")
        done
        
        # Run in parallel with controlled concurrency
        printf '%s\n' "${commands[@]}" | parallel -j "$((THREADS < ${#species_list[@]} ? THREADS : ${#species_list[@]}))" --bar
        
        log_success "Parallel $func_name completed for all species"
    else
        # Fallback to sequential processing
        run_parallel_analysis_fallback "$func_name" "${species_list[@]}"
    fi
}

# Alternative implementation for systems without GNU parallel
run_background_jobs() {
    local func_name=$1
    shift
    local species_list=("$@")
    local max_jobs=${THREADS:-4}
    
    log_info "Running $func_name with background jobs (max: $max_jobs) for ${#species_list[@]} species"
    
    local job_count=0
    local pids=()
    
    for species in "${species_list[@]}"; do
        # Wait if we've reached max jobs
        while [ ${#pids[@]} -ge $max_jobs ]; do
            for i in "${!pids[@]}"; do
                if ! kill -0 "${pids[$i]}" 2>/dev/null; then
                    unset pids[$i]
                fi
            done
            pids=("${pids[@]}")  # Reindex array
            sleep 1
        done
        
        # Start new job
        eval "$func_name $species" &
        pids+=($!)
        log_info "Started $func_name for $species (PID: $!)"
    done
    
    # Wait for all jobs to complete
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
    
    log_success "Background jobs $func_name completed for all species"
}

# Data compression and archiving
compress_outputs() {
    if [ "$compress_outputs" != "true" ]; then
        log_info "Skipping output compression (disabled in config)"
        return 0
    fi
    
    log_info "Compressing large output files..."
    
    # Compress large intermediate files
    find "$OUTDIR/alignments" -name "*.bam" -exec pigz -p "$THREADS" {} \; || true
    find "$OUTDIR/te_counts" -name "*.tsv" -exec pigz -p "$THREADS" {} \; || true
    
    # Create archive of key results
    tar -czf "$OUTDIR/key_results_$(date +%Y%m%d).tar.gz" \
        "$OUTDIR/figures/" \
        "$OUTDIR/reports/" \
        "$OUTDIR/integration_*.rds" \
        "$OUTDIR/orthology/" \
        "$OUTDIR/software_versions.txt"
    
    log_success "Output compression completed"
}

# Performance benchmarking
benchmark_performance() {
    local end_time=$(date +%s)
    local total_time=$((end_time - START_TIME))
    local hours=$((total_time / 3600))
    local minutes=$(((total_time % 3600) / 60))
    local seconds=$((total_time % 60))
    
    log_info "Pipeline completed in ${hours}h ${minutes}m ${seconds}s"
    
    # Generate performance report
    cat > "$OUTDIR/performance_report.txt" << EOF
# Pipeline Performance Report
# Generated on: $(date)

Total Runtime: ${hours}h ${minutes}m ${seconds}s
Species Analyzed: ${#SPECIES_ARRAY[@]}
Threads Used: $THREADS
Memory Allocated: $MEMORY

# Resource Usage Summary
$(cat "$LOG_DIR"/resources_*.log | tail -n +2 | awk -F',' '
BEGIN {max_cpu=0; max_mem=0; total_samples=0}
{
    if($3 > max_cpu) max_cpu = $3
    if($4 > max_mem) max_mem = $4
    total_samples++
}
END {
    print "Peak CPU Usage: " max_cpu "%"
    print "Peak Memory Usage: " max_mem "%"
    print "Monitoring Samples: " total_samples
}')

# Disk Usage
Input Data Size: $(du -sh "$REFDIR" "$FASTQDIR" 2>/dev/null | awk '{sum+=$1} END {print sum "GB"}')
Output Data Size: $(du -sh "$OUTDIR" 2>/dev/null | awk '{print $1}')

EOF
    
    log_success "Performance report generated"
}

# Validation and sanity checks
validate_results() {
    log_info "Validating pipeline results..."
    
    local validation_errors=0
    
    # Check that all expected outputs exist
    local expected_files=(
        "$OUTDIR/software_versions.txt"
        "$OUTDIR/orthology/processed_orthology.rds"
        "$OUTDIR/figures/"
    )
    
    for species in "${SPECIES_ARRAY[@]}"; do
        expected_files+=(
            "$OUTDIR/integration_${species}.rds"
            "$OUTDIR/clusters/advanced_${species}/clustered_data.rds"
            "$OUTDIR/te_counts/harmonized_te_counts_${species}.rds"
        )
    done
    
    for file in "${expected_files[@]}"; do
        if [ ! -e "$file" ]; then
            log_error "Missing expected output: $file"
            ((validation_errors++))
        fi
    done
    
    # Run R-based validation
    Rscript "$SCRIPTSDIR/validate_pipeline_results.R" \
        --pipeline_dir "$OUTDIR" \
        --species_list "${SPECIES_ARRAY[*]}" \
        --output "$OUTDIR/validation_report.txt"
    
    if [ $validation_errors -eq 0 ]; then
        log_success "All validation checks passed"
        return 0
    else
        log_error "$validation_errors validation errors found"
        return 1
    fi
}

# Main pipeline orchestration
main() {
    echo "========================================="
    echo "Enhanced Drosophila TE Analysis Pipeline"
    echo "Version: $PIPELINE_VERSION"
    echo "Started: $(date)"
    echo "========================================="
    echo ""
    
    # Setup
    mkdir -p "$LOG_DIR" "$TEMP_DIR" "$OUTDIR"
    
    # Load configuration
    load_config "$1"
    
    # Load modules if available (for HPC systems)
    if [ -f "load_modules.sh" ]; then
        log_info "Loading HPC modules..."
        source load_modules.sh
    fi
    
    # Check for TE Pipeline environment
    if [ -d "$HOME/te_pipeline" ] && [ -f "$HOME/te_pipeline/activate.sh" ]; then
        log_info "Activating TE Pipeline environment..."
        source "$HOME/te_pipeline/activate.sh"
    fi
    
    # System checks
    check_system_requirements
    record_versions
    
    # Create format detection file
    mkdir -p "$TEMP_DIR"
    touch "$TEMP_DIR/gtf_formats.txt"
    
    echo "Configuration Summary:"
    echo "- Species: ${SPECIES_ARRAY[*]}"
    echo "- Threads: $THREADS"
    echo "- Memory: $MEMORY"
    echo "- Output: $OUTDIR"
    echo ""
    
    # Input validation
    log_info "=== Input Validation ==="
    local validation_failed=false
    for species in "${SPECIES_ARRAY[@]}"; do
        if ! check_files "$species"; then
            validation_failed=true
        fi
    done
    
    if [ "$validation_failed" = true ]; then
        log_error "Input validation failed. Please check your input files."
        exit 1
    fi
    
    # Quality Control Phase
    log_info "=== Quality Control Phase ==="
    run_parallel_analysis "run_fastqc" "${SPECIES_ARRAY[@]}"
    
    # Reference Preparation Phase
    log_info "=== Reference Preparation Phase ==="
    setup_enhanced_orthology
    process_enhanced_orthology
    
    # Alignment Phase
    log_info "=== Alignment Phase ==="
    if command -v parallel &> /dev/null; then
        run_parallel_analysis "build_star_index" "${SPECIES_ARRAY[@]}"
        run_parallel_analysis "run_starsolo" "${SPECIES_ARRAY[@]}"
    else
        log_warning "GNU parallel not available, using background jobs"
        run_background_jobs "build_star_index" "${SPECIES_ARRAY[@]}"
        run_background_jobs "run_starsolo" "${SPECIES_ARRAY[@]}"
    fi
    
    # Quantification Phase
    log_info "=== Quantification Phase ==="
    if command -v parallel &> /dev/null; then
        run_parallel_analysis "run_te_quantification" "${SPECIES_ARRAY[@]}"
    else
        run_background_jobs "run_te_quantification" "${SPECIES_ARRAY[@]}"
    fi
    
    # Single Cell Analysis Phase
    log_info "=== Single Cell Analysis Phase ==="
    if command -v parallel &> /dev/null; then
        run_parallel_analysis "run_cell_qc" "${SPECIES_ARRAY[@]}"
        run_parallel_analysis "run_advanced_clustering" "${SPECIES_ARRAY[@]}"
    else
        run_background_jobs "run_cell_qc" "${SPECIES_ARRAY[@]}"
        run_background_jobs "run_advanced_clustering" "${SPECIES_ARRAY[@]}"
    fi
    
    # Advanced Analysis Phase
    log_info "=== Advanced Analysis Phase ==="
    if command -v parallel &> /dev/null; then
        run_parallel_analysis "run_trajectory_analysis" "${SPECIES_ARRAY[@]}"
        run_parallel_analysis "run_coexpression_analysis" "${SPECIES_ARRAY[@]}"
    else
        run_background_jobs "run_trajectory_analysis" "${SPECIES_ARRAY[@]}"
        run_background_jobs "run_coexpression_analysis" "${SPECIES_ARRAY[@]}"
    fi
    
    # Integration Phase
    log_info "=== Integration Phase ==="
    if command -v parallel &> /dev/null; then
        run_parallel_analysis "integrate_te_clusters_enhanced" "${SPECIES_ARRAY[@]}"
    else
        run_background_jobs "integrate_te_clusters_enhanced" "${SPECIES_ARRAY[@]}"
    fi
    
    # Comparative Analysis Phase
    log_info "=== Comparative Analysis Phase ==="
    run_differential_analysis
    run_evolutionary_analysis
    
    # Visualization Phase
    log_info "=== Visualization Phase ==="
    create_publication_figures
    generate_html_report
    
    # Finalization Phase
    log_info "=== Finalization Phase ==="
    validate_results
    compress_outputs
    benchmark_performance
    
    # Final summary
    echo ""
    echo "========================================="
    echo "Pipeline Completed Successfully!"
    echo "========================================="
    echo ""
    echo "📊 Results Summary:"
    echo "   • Species analyzed: ${#SPECIES_ARRAY[@]}"
    echo "   • Output directory: $OUTDIR"
    echo "   • Key results archive: $OUTDIR/key_results_$(date +%Y%m%d).tar.gz"
    echo ""
    echo "📁 Main Output Directories:"
    echo "   • Quality Control: $OUTDIR/qc/"
    echo "   • Alignments: $OUTDIR/alignments/"
    echo "   • TE Quantification: $OUTDIR/te_counts/"
    echo "   • Cell Clustering: $OUTDIR/clusters/"
    echo "   • Integration Results: $OUTDIR/integration_*.rds"
    echo "   • Comparative Analysis: $OUTDIR/differential_analysis/"
    echo "   • Evolutionary Analysis: $OUTDIR/evolutionary_analysis/"
    echo "   • Publication Figures: $OUTDIR/figures/"
    echo "   • HTML Report: $OUTDIR/reports/pipeline_report.html"
    echo ""
    echo "📈 Performance:"
    local end_time=$(date +%s)
    local total_time=$((end_time - START_TIME))
    local hours=$((total_time / 3600))
    local minutes=$(((total_time % 3600) / 60))
    echo "   • Total runtime: ${hours}h ${minutes}m"
    echo "   • Performance report: $OUTDIR/performance_report.txt"
    echo ""
    echo "✅ Pipeline completed at: $(date)"
    echo ""
    echo "🔬 Ready for publication! Check the HTML report for detailed results."
    echo "========================================="
}

# Helper function for usage information
show_usage() {
    cat << EOF
Enhanced Drosophila TE Analysis Pipeline v$PIPELINE_VERSION

USAGE:
    $0 [CONFIG_FILE]

ARGUMENTS:
    CONFIG_FILE    Path to pipeline configuration file (optional)
                  Default: $DEFAULT_CONFIG

CONFIGURATION:
    The pipeline uses a YAML-style configuration file to specify parameters.
    If no config file is provided, a default one will be created.

REQUIREMENTS:
    • STAR aligner
    • samtools
    • R with required packages (Seurat, SAKE, etc.)
    • Python3 with required packages
    • telescope
    • featureCounts
    • parallel
    • pigz (optional, for compression)

INPUT STRUCTURE:
    references/
    ├── dmel_genome.fasta
    ├── dmel_genes.harmonized.gtf
    ├── dmelTEs.harmonized.gtf
    ├── dyak_genome.fasta
    ├── dyak_genes.harmonized.gtf
    ├── dyakTEs.harmonized.gtf
    ├── dana_genome.fasta
    ├── dana_genes.harmonized.gtf
    └── danaTEs.harmonized.gtf
    
    fastqs/
    ├── Dmel/
    │   ├── sample1_R1.fastq.gz
    │   └── sample1_R2.fastq.gz
    ├── Dyak/
    └── Dana/

OUTPUT STRUCTURE:
    results/
    ├── qc/                    # Quality control reports
    ├── alignments/            # STAR alignments
    ├── te_counts/            # TE quantification
    ├── clusters/             # Cell clustering results
    ├── integration_*.rds     # Integrated TE-cluster data
    ├── differential_analysis/ # Cross-species comparisons
    ├── evolutionary_analysis/ # Evolutionary insights
    ├── figures/              # Publication-ready figures
    ├── reports/              # HTML reports
    └── software_versions.txt # Reproducibility info

EXAMPLES:
    # Run with default configuration
    $0
    
    # Run with custom configuration
    $0 my_config.yaml
    
    # Create a default config file
    mkdir -p config && $0 config/my_config.yaml

For more information, see the generated HTML report after pipeline completion.
EOF
}

# Script entry point
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            -v|--version)
                echo "Enhanced Drosophila TE Analysis Pipeline v$PIPELINE_VERSION"
                exit 0
                ;;
            *)
                DEFAULT_CONFIG="$1"
                shift
                ;;
        esac
    done
    
    # Run the main pipeline
    main "$DEFAULT_CONFIG"
fi\t'; then
                log_error "Invalid GTF format in $gtf_file"
                ((errors++))
            else
                local feature_count=$(wc -l < "$gtf_file")
                log_info "GTF file $gtf_file contains $feature_count features (format: $gtf_format)"
                
                # Store format information for later use
                if [[ "$gtf_file" == *"genes"* ]] || [[ "$gtf_file" == *".gtf" ]] && [[ "$gtf_file" != *"TE"* ]]; then
                    echo "$species:genes:$gtf_format" >> "$TEMP_DIR/gtf_formats.txt"
                elif [[ "$gtf_file" == *"TE"* ]] || [[ "$gtf_file" == *"repeat"* ]]; then
                    echo "$species:tes:$gtf_format" >> "$TEMP_DIR/gtf_formats.txt"
                fi
            fi
        fi
    done
    
    # Check FASTQ directory and files
    if [ ! -d "$fastq_dir" ]; then
        log_error "Missing FASTQ directory: $fastq_dir"
        ((errors++))
    else
        local fastq_files=($(find "$fastq_dir" -name "*.fastq.gz" | sort))
        if [ ${#fastq_files[@]} -lt 2 ]; then
            log_error "Expected at least 2 FASTQ files in $fastq_dir, found ${#fastq_files[@]}"
            ((errors++))
        else
            log_info "Found ${#fastq_files[@]} FASTQ files for $species"
            
            # Validate FASTQ format and check for paired-end
            for fq in "${fastq_files[@]:0:2}"; do
                if ! zcat "$fq" 2>/dev/null | head -4 | grep -q "^@"; then
                    log_error "Invalid FASTQ format in $fq"
                    ((errors++))
                fi
            done
        fi
    fi
    
    if [ $errors -eq 0 ]; then
        log_success "All files validated for $species"
        return 0
    else
        log_error "$errors validation errors found for $species"
        return 1
    fi
}

# Resource monitoring function
monitor_resources() {
    local pid=$1
    local process_name=$2
    local logfile="$LOG_DIR/resources_${process_name}.log"
    
    echo "Time,PID,CPU%,MEM%,VSZ,RSS" > "$logfile"
    
    while kill -0 $pid 2>/dev/null; do
        ps -p $pid -o pid,pcpu,pmem,vsz,rss --no-headers >> "$logfile" 2>/dev/null || break
        sleep 30
    done &
    
    echo $!  # Return monitoring PID
}

# Enhanced quality control with FastQC
run_fastqc() {
    local species=$1
    log_info "Running FastQC quality control for $species..."
    
    local fastq_dir="$FASTQDIR/$species"
    local qc_dir="$OUTDIR/qc/fastqc_${species}"
    
    mkdir -p "$qc_dir"
    
    # Run FastQC on all FASTQ files
    fastqc -t "$THREADS" -o "$qc_dir" "$fastq_dir"/*.fastq.gz
    
    # Generate summary report
    python3 "$SCRIPTSDIR/summarize_fastqc.py" \
        --input_dir "$qc_dir" \
        --output "$qc_dir/fastqc_summary.json"
    
    log_success "FastQC completed for $species"
}

# Version tracking
record_versions() {
    local version_file="$OUTDIR/software_versions.txt"
    
    log_info "Recording software versions..."
    
    cat > "$version_file" << EOF
# Software Versions Used in Analysis
# Generated on: $(date)
# Pipeline Version: $PIPELINE_VERSION

EOF
    
    # Record versions of key tools
    echo "STAR: $(STAR --version 2>&1 | head -1)" >> "$version_file"
    echo "samtools: $(samtools --version | head -1)" >> "$version_file"
    echo "R: $(Rscript --version 2>&1 | head -1)" >> "$version_file"
    echo "Python: $(python3 --version)" >> "$version_file"
    
    # Record R package versions
    Rscript -e "
    packages <- c('Seurat', 'dplyr', 'ggplot2', 'SingleCellExperiment', 'scater', 'SAKE', 'DESeq2', 'monocle3')
    for(pkg in packages) {
        if(require(pkg, quietly=TRUE, character.only=TRUE)) {
            cat(paste0(pkg, ': ', packageVersion(pkg), '\n'))
        }
    }
    " >> "$version_file"
    
    log_success "Software versions recorded"
}

# Enhanced STAR index building with validation
build_star_index() {
    local species=$1
    log_info "Building STAR index for $species..."
    
    local genome_file="$REFDIR/${species,,}_genome.fasta"
    local genes_file="$REFDIR/${species,,}_genes.harmonized.gtf"
    local index_dir="$OUTDIR/star_indices/star_index_${species}"
    
    mkdir -p "$index_dir"
    
    # Calculate optimal parameters based on genome size
    local genome_size=$(stat -c%s "$genome_file")
    local genomeSAindexNbases=$((14))  # Default for most genomes
    
    if [ $genome_size -lt 2000000000 ]; then  # < 2GB
        genomeSAindexNbases=13
    fi
    
    # Run STAR with resource monitoring
    STAR --runMode genomeGenerate \
         --genomeDir "$index_dir" \
         --genomeFastaFiles "$genome_file" \
         --sjdbGTFfile "$genes_file" \
         --sjdbOverhang 100 \
         --runThreadN "$THREADS" \
         --limitGenomeGenerateRAM $(($(echo $MEMORY | sed 's/G//') * 1000000000)) \
         --genomeSAindexNbases $genomeSAindexNbases \
         --outTmpDir "$TEMP_DIR/star_tmp_${species}" &
    
    local star_pid=$!
    local monitor_pid=$(monitor_resources $star_pid "star_index_${species}")
    
    wait $star_pid
    kill $monitor_pid 2>/dev/null || true
    
    # Validate index
    if [ ! -f "$index_dir/SA" ]; then
        log_error "STAR index building failed for $species"
        return 1
    fi
    
    log_success "STAR index built successfully for $species"
}

# Enhanced STARsolo with comprehensive parameters
run_starsolo() {
    local species=$1
    log_info "Running STARsolo alignment for $species..."
    
    local fastq_dir="$FASTQDIR/$species"
    local fastq_files=($(find "$fastq_dir" -name "*.fastq.gz" | sort))
    local fq1="${fastq_files[0]}"
    local fq2="${fastq_files[1]}"
    local index_dir="$OUTDIR/star_indices/star_index_${species}"
    local output_dir="$OUTDIR/alignments/starsolo_${species}"
    
    mkdir -p "$output_dir"
    
    # Enhanced STARsolo parameters
    STAR --genomeDir "$index_dir" \
         --readFilesIn "$fq1" "$fq2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$output_dir/" \
         --soloType CB_UMI_Simple \
         --soloCBwhitelist None \
         --soloCBstart 1 --soloCBlen 16 \
         --soloUMIstart 17 --soloUMIlen 12 \
         --soloStrand Forward \
         --soloFeatures Gene GeneFull Velocyto \
         --soloOutFileNames Solo.out/ features.tsv barcodes.tsv matrix.mtx \
         --soloMultiMappers EM \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --runThreadN "$THREADS" \
         --limitBAMsortRAM $(($(echo $MEMORY | sed 's/G//') * 1000000000)) \
         --outTmpDir "$TEMP_DIR/starsolo_tmp_${species}" &
    
    local star_pid=$!
    local monitor_pid=$(monitor_resources $star_pid "starsolo_${species}")
    
    wait $star_pid
    kill $monitor_pid 2>/dev/null || true
    
    # Copy and index the sorted BAM file
    cp "$output_dir/Aligned.sortedByCoord.out.bam" "$OUTDIR/alignments/starsolo_${species}.bam"
    samtools index "$OUTDIR/alignments/starsolo_${species}.bam"
    
    # Generate alignment statistics
    collect_alignment_stats "$species"
    
    log_success "STARsolo alignment completed for $species"
}

# Comprehensive alignment statistics
collect_alignment_stats() {
    local species=$1
    log_info "Collecting alignment statistics for $species..."
    
    local bam_file="$OUTDIR/alignments/starsolo_${species}.bam"
    local stats_dir="$OUTDIR/qc/alignment_stats"
    mkdir -p "$stats_dir"
    
    # Basic alignment statistics
    samtools flagstat "$bam_file" > "$stats_dir/${species}_flagstat.txt"
    samtools idxstats "$bam_file" > "$stats_dir/${species}_idxstats.txt"
    samtools stats "$bam_file" > "$stats_dir/${species}_stats.txt"
    
    # STARsolo statistics
    local starsolo_dir="$OUTDIR/alignments/starsolo_${species}"
    if [ -f "$starsolo_dir/Log.final.out" ]; then
        cp "$starsolo_dir/Log.final.out" "$stats_dir/${species}_star_final.log"
    fi
    
    # Generate summary statistics
    python3 "$SCRIPTSDIR/summarize_alignment_stats.py" \
        --species "$species" \
        --stats_dir "$stats_dir" \
        --output "$stats_dir/${species}_alignment_summary.json"
    
    log_success "Alignment statistics collected for $species"
}

# Enhanced TE quantification with multiple methods
run_te_quantification() {
    local species=$1
    log_info "Running comprehensive TE quantification for $species..."
    
    local bam_file="$OUTDIR/alignments/starsolo_${species}.bam"
    local te_file="$REFDIR/${species,,}TEs.harmonized.gtf"
    local output_dir="$OUTDIR/te_counts"
    
    mkdir -p "$output_dir"
    
    # Sort BAM by name for Telescope
    local name_sorted_bam="$OUTDIR/alignments/starsolo_${species}.sortedByName.bam"
    samtools sort -n -@ "$THREADS" -o "$name_sorted_bam" "$bam_file"
    
    # Run Telescope
    log_info "Running Telescope for $species..."
    telescope assign \
        --barcode_tag CB \
        --umi_tag UB \
        --stranded_library \
        --theta_prior 200000 \
        --max_iter 200 \
        --updated_sam \
        "$name_sorted_bam" \
        "$te_file" \
        > "$output_dir/telescope_${species}_report.tsv" &
    
    local telescope_pid=$!
    
    # Run featureCounts for comparison
    log_info "Running featureCounts for $species..."
    featureCounts \
        -T "$THREADS" \
        -p \
        -B \
        -C \
        -s 1 \
        -a "$te_file" \
        -o "$output_dir/featurecounts_${species}_te.txt" \
        "$bam_file" &
    
    local fc_pid=$!
    
    # Wait for both methods
    wait $telescope_pid $fc_pid
    
    # Process and harmonize results
    Rscript "$SCRIPTSDIR/harmonize_te_counts.R" \
        --telescope "$output_dir/telescope_${species}_report.tsv" \
        --featurecounts "$output_dir/featurecounts_${species}_te.txt" \
        --species "$species" \
        --output "$output_dir/harmonized_te_counts_${species}.rds"
    
    log_success "TE quantification completed for $species"
}

# Advanced cell filtering and quality control
run_cell_qc() {
    local species=$1
    log_info "Running advanced cell quality control for $species..."
    
    local starsolo_dir="$OUTDIR/alignments/starsolo_${species}"
    local qc_dir="$OUTDIR/qc/cell_qc_${species}"
    
    mkdir -p "$qc_dir"
    
    Rscript "$SCRIPTSDIR/advanced_cell_qc.R" \
        --input "$starsolo_dir/Solo.out/Gene/filtered/" \
        --species "$species" \
        --min_cells "$MIN_CELLS" \
        --min_genes "$MIN_GENES" \
        --max_mito_percent "$MAX_MITO_PERCENT" \
        --output "$qc_dir/" \
        --seed "$SEED"
    
    log_success "Cell quality control completed for $species"
}

# Enhanced clustering with multiple algorithms
run_advanced_clustering() {
    local species=$1
    log_info "Running advanced clustering analysis for $species..."
    
    local qc_dir="$OUTDIR/qc/cell_qc_${species}"
    local cluster_dir="$OUTDIR/clusters/advanced_${species}"
    
    mkdir -p "$cluster_dir"
    
    # Run multiple clustering algorithms
    Rscript "$SCRIPTSDIR/advanced_clustering.R" \
        --input "$qc_dir/filtered_data.rds" \
        --species "$species" \
        --resolution "$RESOLUTION" \
        --algorithms "leiden,louvain,hierarchical" \
        --output "$cluster_dir/" \
        --seed "$SEED"
    
    log_success "Advanced clustering completed for $species"
}

# Trajectory analysis
run_trajectory_analysis() {
    local species=$1
    
    if [ "$run_trajectory" != "true" ]; then
        log_info "Skipping trajectory analysis for $species (disabled in config)"
        return 0
    fi
    
    log_info "Running trajectory analysis for $species..."
    
    local cluster_dir="$OUTDIR/clusters/advanced_${species}"
    local trajectory_dir="$OUTDIR/trajectory/${species}"
    
    mkdir -p "$trajectory_dir"
    
    Rscript "$SCRIPTSDIR/trajectory_analysis.R" \
        --input "$cluster_dir/clustered_data.rds" \
        --species "$species" \
        --output "$trajectory_dir/" \
        --seed "$SEED"
    
    log_success "Trajectory analysis completed for $species"
}

# Co-expression network analysis
run_coexpression_analysis() {
    local species=$1
    
    if [ "$run_coexpression" != "true" ]; then
        log_info "Skipping co-expression analysis for $species (disabled in config)"
        return 0
    fi
    
    log_info "Running co-expression network analysis for $species..."
    
    local te_counts="$OUTDIR/te_counts/harmonized_te_counts_${species}.rds"
    local gene_counts="$OUTDIR/alignments/starsolo_${species}/Solo.out/Gene/filtered/"
    local coexp_dir="$OUTDIR/coexpression/${species}"
    
    mkdir -p "$coexp_dir"
    
    Rscript "$SCRIPTSDIR/coexpression_analysis.R" \
        --te_counts "$te_counts" \
        --gene_counts "$gene_counts" \
        --species "$species" \
        --output "$coexp_dir/" \
        --seed "$SEED"
    
    log_success "Co-expression analysis completed for $species"
}

# Enhanced orthology analysis with phylogenetic context
setup_enhanced_orthology() {
    log_info "Setting up enhanced orthology analysis..."
    
    local ortho_dir="$REFDIR/orthology"
    mkdir -p "$ortho_dir"
    cd "$ortho_dir"
    
    # Download FlyBase orthology data
    if [ ! -f "dmel_orthologs_in_drosophila_species_fb_2023_05.tsv" ]; then
        log_info "Downloading FlyBase orthology data..."
        wget -O dmel_orthologs_in_drosophila_species_fb_2023_05.tsv.gz \
            "$flybase_orthology_url"
        gunzip dmel_orthologs_in_drosophila_species_fb_2023_05.tsv.gz
    fi
    
    # Download Dfam data for TE classification
    if [ ! -f "Dfam_3.7_families.tsv" ]; then
        log_info "Downloading Dfam TE classification data..."
        wget -O Dfam_3.7_families.tsv.gz "$dfam_url"
        gunzip Dfam_3.7_families.tsv.gz
    fi
    
    # Download phylogenetic tree
    if [ ! -f "drosophila_species_tree.nwk" ]; then
        log_info "Creating phylogenetic tree for Drosophila species..."
        cat > drosophila_species_tree.nwk << 'EOF'
(Dmel:0.1,(Dyak:0.05,Dana:0.05):0.05);
EOF
    fi
    
    cd - > /dev/null
    log_success "Enhanced orthology setup completed"
}

# Comprehensive orthology processing with evolutionary analysis
process_enhanced_orthology() {
    log_info "Processing enhanced orthology data with evolutionary analysis..."
    
    Rscript "$SCRIPTSDIR/enhanced_orthology_analysis.R" \
        --orthology_file "$REFDIR/orthology/dmel_orthologs_in_drosophila_species_fb_2023_05.tsv" \
        --dfam_file "$REFDIR/orthology/Dfam_3.7_families.tsv" \
        --tree_file "$REFDIR/orthology/drosophila_species_tree.nwk" \
        --species_list "${SPECIES_ARRAY[*]}" \
        --output "$OUTDIR/orthology/"
    
    log_success "Enhanced orthology processing completed"
}

# Comprehensive differential analysis
run_differential_analysis() {
    if [ "$run_differential" != "true" ]; then
        log_info "Skipping differential analysis (disabled in config)"
        return 0
    fi
    
    log_info "Running comprehensive differential analysis..."
    
    local diff_dir="$OUTDIR/differential_analysis"
    mkdir -p "$diff_dir"
    
    # Collect all integration files
    local integration_files=()
    for species in "${SPECIES_ARRAY[@]}"; do
        integration_files+=("$OUTDIR/integration_${species}.rds")
    done
    
    Rscript "$SCRIPTSDIR/comprehensive_differential_analysis.R" \
        --integration_files "${integration_files[@]}" \
        --orthology_data "$OUTDIR/orthology/processed_orthology.rds" \
        --output "$diff_dir/" \
        --species_list "${SPECIES_ARRAY[*]}" \
        --seed "$SEED"
    
    log_success "Differential analysis completed"
}

# Evolutionary analysis
run_evolutionary_analysis() {
    if [ "$run_evolutionary" != "true" ]; then
        log_info "Skipping evolutionary analysis (disabled in config)"
        return 0
    fi
    
    log_info "Running evolutionary analysis..."
    
    local evo_dir="$OUTDIR/evolutionary_analysis"
    mkdir -p "$evo_dir"
    
    Rscript "$SCRIPTSDIR/evolutionary_analysis.R" \
        --integration_dir "$OUTDIR/" \
        --orthology_data "$OUTDIR/orthology/processed_orthology.rds" \
        --tree_file "$REFDIR/orthology/drosophila_species_tree.nwk" \
        --output "$evo_dir/" \
        --species_list "${SPECIES_ARRAY[*]}"
    
    log_success "Evolutionary analysis completed"
}

# Enhanced integration of TE counts with clusters
integrate_te_clusters_enhanced() {
    local species=$1
    log_info "Running enhanced TE-cluster integration for $species..."
    
    local te_counts="$OUTDIR/te_counts/harmonized_te_counts_${species}.rds"
    local cluster_data="$OUTDIR/clusters/advanced_${species}/clustered_data.rds"
    local cell_qc="$OUTDIR/qc/cell_qc_${species}/cell_metadata.rds"
    local output_file="$OUTDIR/integration_${species}.rds"
    
    Rscript "$SCRIPTSDIR/enhanced_te_cluster_integration.R" \
        --te_counts "$te_counts" \
        --cluster_data "$cluster_data" \
        --cell_qc "$cell_qc" \
        --species "$species" \
        --output "$output_file" \
        --min_te_count "$MIN_TE_COUNT"
    
    log_success "Enhanced TE-cluster integration completed for $species"
}

# Publication-ready figure generation
create_publication_figures() {
    log_info "Creating publication-ready figures..."
    
    local figures_dir="$OUTDIR/figures"
    mkdir -p "$figures_dir"
    
    # Collect all integration files
    local integration_files=()
    for species in "${SPECIES_ARRAY[@]}"; do
        if [ -f "$OUTDIR/integration_${species}.rds" ]; then
            integration_files+=("$OUTDIR/integration_${species}.rds")
        fi
    done
    
    Rscript "$SCRIPTSDIR/create_publication_figures.R" \
        --integration_files "${integration_files[@]}" \
        --orthology_data "$OUTDIR/orthology/processed_orthology.rds" \
        --differential_data "$OUTDIR/differential_analysis/" \
        --evolutionary_data "$OUTDIR/evolutionary_analysis/" \
        --output "$figures_dir/" \
        --species_list "${SPECIES_ARRAY[*]}"
    
    log_success "Publication-ready figures created"
}

# Comprehensive HTML report generation
generate_html_report() {
    if [ "$create_html_report" != "true" ]; then
        log_info "Skipping HTML report generation (disabled in config)"
        return 0
    fi
    
    log_info "Generating comprehensive HTML report..."
    
    local report_dir="$OUTDIR/reports"
    mkdir -p "$report_dir"
    
    Rscript "$SCRIPTSDIR/generate_html_report.R" \
        --pipeline_dir "$OUTDIR" \
        --config_file "$DEFAULT_CONFIG" \
        --versions_file "$OUTDIR/software_versions.txt" \
        --output "$report_dir/pipeline_report.html" \
        --species_list "${SPECIES_ARRAY[*]}"
    
    log_success "HTML report generated at $report_dir/pipeline_report.html"
}

# Parallel processing wrapper
run_parallel_analysis() {
    local func_name=$1
    shift
    local species_list=("$@")
    
    log_info "Running $func_name in parallel for ${#species_list[@]} species"
    
    # Create array of commands
    local commands=()
    for species in "${species_list[@]}"; do
        commands+=("$func_name $species")
    done
    
    # Run in parallel with controlled concurrency
    printf '%s\n' "${commands[@]}" | parallel -j "$THREADS" --bar
    
    log_success "Parallel $func_name completed for all species"
}

# Data compression and archiving
compress_outputs() {
    if [ "$compress_outputs" != "true" ]; then
        log_info "Skipping output compression (disabled in config)"
        return 0
    fi
    
    log_info "Compressing large output files..."
    
    # Compress large intermediate files
    find "$OUTDIR/alignments" -name "*.bam" -exec pigz -p "$THREADS" {} \; || true
    find "$OUTDIR/te_counts" -name "*.tsv" -exec pigz -p "$THREADS" {} \; || true
    
    # Create archive of key results
    tar -czf "$OUTDIR/key_results_$(date +%Y%m%d).tar.gz" \
        "$OUTDIR/figures/" \
        "$OUTDIR/reports/" \
        "$OUTDIR/integration_*.rds" \
        "$OUTDIR/orthology/" \
        "$OUTDIR/software_versions.txt"
    
    log_success "Output compression completed"
}

# Performance benchmarking
benchmark_performance() {
    local end_time=$(date +%s)
    local total_time=$((end_time - START_TIME))
    local hours=$((total_time / 3600))
    local minutes=$(((total_time % 3600) / 60))
    local seconds=$((total_time % 60))
    
    log_info "Pipeline completed in ${hours}h ${minutes}m ${seconds}s"
    
    # Generate performance report
    cat > "$OUTDIR/performance_report.txt" << EOF
# Pipeline Performance Report
# Generated on: $(date)

Total Runtime: ${hours}h ${minutes}m ${seconds}s
Species Analyzed: ${#SPECIES_ARRAY[@]}
Threads Used: $THREADS
Memory Allocated: $MEMORY

# Resource Usage Summary
$(cat "$LOG_DIR"/resources_*.log | tail -n +2 | awk -F',' '
BEGIN {max_cpu=0; max_mem=0; total_samples=0}
{
    if($3 > max_cpu) max_cpu = $3
    if($4 > max_mem) max_mem = $4
    total_samples++
}
END {
    print "Peak CPU Usage: " max_cpu "%"
    print "Peak Memory Usage: " max_mem "%"
    print "Monitoring Samples: " total_samples
}')

# Disk Usage
Input Data Size: $(du -sh "$REFDIR" "$FASTQDIR" 2>/dev/null | awk '{sum+=$1} END {print sum "GB"}')
Output Data Size: $(du -sh "$OUTDIR" 2>/dev/null | awk '{print $1}')

EOF
    
    log_success "Performance report generated"
}

# Validation and sanity checks
validate_results() {
    log_info "Validating pipeline results..."
    
    local validation_errors=0
    
    # Check that all expected outputs exist
    local expected_files=(
        "$OUTDIR/software_versions.txt"
        "$OUTDIR/orthology/processed_orthology.rds"
        "$OUTDIR/figures/"
    )
    
    for species in "${SPECIES_ARRAY[@]}"; do
        expected_files+=(
            "$OUTDIR/integration_${species}.rds"
            "$OUTDIR/clusters/advanced_${species}/clustered_data.rds"
            "$OUTDIR/te_counts/harmonized_te_counts_${species}.rds"
        )
    done
    
    for file in "${expected_files[@]}"; do
        if [ ! -e "$file" ]; then
            log_error "Missing expected output: $file"
            ((validation_errors++))
        fi
    done
    
    # Run R-based validation
    Rscript "$SCRIPTSDIR/validate_pipeline_results.R" \
        --pipeline_dir "$OUTDIR" \
        --species_list "${SPECIES_ARRAY[*]}" \
        --output "$OUTDIR/validation_report.txt"
    
    if [ $validation_errors -eq 0 ]; then
        log_success "All validation checks passed"
        return 0
    else
        log_error "$validation_errors validation errors found"
        return 1
    fi
}

# Main pipeline orchestration
main() {
    echo "========================================="
    echo "Enhanced Drosophila TE Analysis Pipeline"
    echo "Version: $PIPELINE_VERSION"
    echo "Started: $(date)"
    echo "========================================="
    echo ""
    
    # Setup
    mkdir -p "$LOG_DIR" "$TEMP_DIR" "$OUTDIR"
    
    # Load configuration
    load_config "$1"
    
    # System checks
    check_system_requirements
    record_versions
    
    echo "Configuration Summary:"
    echo "- Species: ${SPECIES_ARRAY[*]}"
    echo "- Threads: $THREADS"
    echo "- Memory: $MEMORY"
    echo "- Output: $OUTDIR"
    echo ""
    
    # Input validation
    log_info "=== Input Validation ==="
    local validation_failed=false
    for species in "${SPECIES_ARRAY[@]}"; do
        if ! check_files "$species"; then
            validation_failed=true
        fi
    done
    
    if [ "$validation_failed" = true ]; then
        log_error "Input validation failed. Please check your input files."
        exit 1
    fi
    
    # Quality Control Phase
    log_info "=== Quality Control Phase ==="
    run_parallel_analysis "run_fastqc" "${SPECIES_ARRAY[@]}"
    
    # Reference Preparation Phase
    log_info "=== Reference Preparation Phase ==="
    setup_enhanced_orthology
    process_enhanced_orthology
    
    # Alignment Phase
    log_info "=== Alignment Phase ==="
    run_parallel_analysis "build_star_index" "${SPECIES_ARRAY[@]}"
    run_parallel_analysis "run_starsolo" "${SPECIES_ARRAY[@]}"
    
    # Quantification Phase
    log_info "=== Quantification Phase ==="
    run_parallel_analysis "run_te_quantification" "${SPECIES_ARRAY[@]}"
    
    # Single Cell Analysis Phase
    log_info "=== Single Cell Analysis Phase ==="
    run_parallel_analysis "run_cell_qc" "${SPECIES_ARRAY[@]}"
    run_parallel_analysis "run_advanced_clustering" "${SPECIES_ARRAY[@]}"
    
    # Advanced Analysis Phase
    log_info "=== Advanced Analysis Phase ==="
    run_parallel_analysis "run_trajectory_analysis" "${SPECIES_ARRAY[@]}"
    run_parallel_analysis "run_coexpression_analysis" "${SPECIES_ARRAY[@]}"
    
    # Integration Phase
    log_info "=== Integration Phase ==="
    run_parallel_analysis "integrate_te_clusters_enhanced" "${SPECIES_ARRAY[@]}"
    
    # Comparative Analysis Phase
    log_info "=== Comparative Analysis Phase ==="
    run_differential_analysis
    run_evolutionary_analysis
    
    # Visualization Phase
    log_info "=== Visualization Phase ==="
    create_publication_figures
    generate_html_report
    
    # Finalization Phase
    log_info "=== Finalization Phase ==="
    validate_results
    compress_outputs
    benchmark_performance
    
    # Final summary
    echo ""
    echo "========================================="
    echo "Pipeline Completed Successfully!"
    echo "========================================="
    echo ""
    echo "📊 Results Summary:"
    echo "   • Species analyzed: ${#SPECIES_ARRAY[@]}"
    echo "   • Output directory: $OUTDIR"
    echo "   • Key results archive: $OUTDIR/key_results_$(date +%Y%m%d).tar.gz"
    echo ""
    echo "📁 Main Output Directories:"
    echo "   • Quality Control: $OUTDIR/qc/"
    echo "   • Alignments: $OUTDIR/alignments/"
    echo "   • TE Quantification: $OUTDIR/te_counts/"
    echo "   • Cell Clustering: $OUTDIR/clusters/"
    echo "   • Integration Results: $OUTDIR/integration_*.rds"
    echo "   • Comparative Analysis: $OUTDIR/differential_analysis/"
    echo "   • Evolutionary Analysis: $OUTDIR/evolutionary_analysis/"
    echo "   • Publication Figures: $OUTDIR/figures/"
    echo "   • HTML Report: $OUTDIR/reports/pipeline_report.html"
    echo ""
    echo "📈 Performance:"
    local end_time=$(date +%s)
    local total_time=$((end_time - START_TIME))
    local hours=$((total_time / 3600))
    local minutes=$(((total_time % 3600) / 60))
    echo "   • Total runtime: ${hours}h ${minutes}m"
    echo "   • Performance report: $OUTDIR/performance_report.txt"
    echo ""
    echo "✅ Pipeline completed at: $(date)"
    echo ""
    echo "🔬 Ready for publication! Check the HTML report for detailed results."
    echo "========================================="
}

# Helper function for usage information
show_usage() {
    cat << EOF
Enhanced Drosophila TE Analysis Pipeline v$PIPELINE_VERSION

USAGE:
    $0 [CONFIG_FILE]

ARGUMENTS:
    CONFIG_FILE    Path to pipeline configuration file (optional)
                  Default: $DEFAULT_CONFIG

CONFIGURATION:
    The pipeline uses a YAML-style configuration file to specify parameters.
    If no config file is provided, a default one will be created.

REQUIREMENTS:
    • STAR aligner
    • samtools
    • R with required packages (Seurat, SAKE, etc.)
    • Python3 with required packages
    • telescope
    • featureCounts
    • parallel
    • pigz (optional, for compression)

INPUT STRUCTURE:
    references/
    ├── dmel_genome.fasta
    ├── dmel_genes.harmonized.gtf
    ├── dmelTEs.harmonized.gtf
    ├── dyak_genome.fasta
    ├── dyak_genes.harmonized.gtf
    ├── dyakTEs.harmonized.gtf
    ├── dana_genome.fasta
    ├── dana_genes.harmonized.gtf
    └── danaTEs.harmonized.gtf
    
    fastqs/
    ├── Dmel/
    │   ├── sample1_R1.fastq.gz
    │   └── sample1_R2.fastq.gz
    ├── Dyak/
    └── Dana/

OUTPUT STRUCTURE:
    results/
    ├── qc/                    # Quality control reports
    ├── alignments/            # STAR alignments
    ├── te_counts/            # TE quantification
    ├── clusters/             # Cell clustering results
    ├── integration_*.rds     # Integrated TE-cluster data
    ├── differential_analysis/ # Cross-species comparisons
    ├── evolutionary_analysis/ # Evolutionary insights
    ├── figures/              # Publication-ready figures
    ├── reports/              # HTML reports
    └── software_versions.txt # Reproducibility info

EXAMPLES:
    # Run with default configuration
    $0
    
    # Run with custom configuration
    $0 my_config.yaml
    
    # Create a default config file
    mkdir -p config && $0 config/my_config.yaml

For more information, see the generated HTML report after pipeline completion.
EOF
}

# Script entry point
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            -v|--version)
                echo "Enhanced Drosophila TE Analysis Pipeline v$PIPELINE_VERSION"
                exit 0
                ;;
            *)
                DEFAULT_CONFIG="$1"
                shift
                ;;
        esac
    done
    
    # Run the main pipeline
    main "$DEFAULT_CONFIG"
fi
