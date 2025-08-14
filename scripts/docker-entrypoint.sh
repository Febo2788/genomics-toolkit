#!/bin/bash
set -e

# GenomicsToolkit Docker Entrypoint Script
# This script handles initialization and execution of the GenomicsToolkit in Docker

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')] $1${NC}"
}

error() {
    echo -e "${RED}[ERROR] $1${NC}" >&2
}

warn() {
    echo -e "${YELLOW}[WARN] $1${NC}"
}

success() {
    echo -e "${GREEN}[SUCCESS] $1${NC}"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to initialize configuration
init_config() {
    log "Initializing GenomicsToolkit configuration..."
    
    if [ ! -f "$GENOMICS_TOOLKIT_CONFIG/config.yaml" ]; then
        log "Creating default configuration file..."
        mkdir -p "$GENOMICS_TOOLKIT_CONFIG"
        
        # Generate default config
        genomics-toolkit config --list-all > /dev/null 2>&1 || {
            warn "Could not generate default config, creating minimal config"
            cat > "$GENOMICS_TOOLKIT_CONFIG/config.yaml" << EOF
logging:
  level: INFO
  format: '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

processing:
  max_workers: 4
  chunk_size: 1000
  memory_limit: 8GB

quality_control:
  min_sequence_length: 50
  max_n_content: 10
  min_quality_score: 20

databases:
  cache_dir: /app/data/cache
  request_delay: 0.5
EOF
        }
        success "Configuration initialized"
    else
        log "Configuration already exists"
    fi
}

# Function to set up directories
setup_directories() {
    log "Setting up directories..."
    
    # Create required directories
    mkdir -p "$GENOMICS_TOOLKIT_DATA"
    mkdir -p "$GENOMICS_TOOLKIT_OUTPUT"
    mkdir -p "$GENOMICS_TOOLKIT_CONFIG"
    mkdir -p /app/data/cache
    mkdir -p /app/logs
    
    # Ensure correct permissions
    if [ "$USER" = "genomics" ]; then
        chown -R genomics:genomics "$GENOMICS_TOOLKIT_OUTPUT" 2>/dev/null || true
        chown -R genomics:genomics "$GENOMICS_TOOLKIT_CONFIG" 2>/dev/null || true
        chown -R genomics:genomics /app/data/cache 2>/dev/null || true
        chown -R genomics:genomics /app/logs 2>/dev/null || true
    fi
    
    success "Directories set up"
}

# Function to check system resources
check_resources() {
    log "Checking system resources..."
    
    # Check available memory
    TOTAL_MEM=$(cat /proc/meminfo | grep MemTotal | awk '{print $2}')
    AVAILABLE_MEM=$(cat /proc/meminfo | grep MemAvailable | awk '{print $2}')
    MEM_GB=$((TOTAL_MEM / 1024 / 1024))
    AVAIL_GB=$((AVAILABLE_MEM / 1024 / 1024))
    
    log "Total memory: ${MEM_GB}GB, Available: ${AVAIL_GB}GB"
    
    if [ $AVAIL_GB -lt 2 ]; then
        warn "Low available memory (${AVAIL_GB}GB). Some operations may be slow."
    fi
    
    # Check CPU count
    CPU_COUNT=$(nproc)
    log "CPU cores available: $CPU_COUNT"
    
    # Check disk space
    DISK_AVAIL=$(df /app/output | tail -1 | awk '{print $4}')
    DISK_GB=$((DISK_AVAIL / 1024 / 1024))
    log "Available disk space: ${DISK_GB}GB"
    
    if [ $DISK_GB -lt 1 ]; then
        warn "Low disk space (${DISK_GB}GB). Ensure sufficient space for output files."
    fi
}

# Function to validate installation
validate_installation() {
    log "Validating GenomicsToolkit installation..."
    
    # Check if genomics-toolkit command is available
    if ! command_exists genomics-toolkit; then
        error "genomics-toolkit command not found!"
        exit 1
    fi
    
    # Check if Python modules can be imported
    python -c "
import genomics_toolkit
from genomics_toolkit import SequenceAnalyzer, VariantCaller, Visualizer
print('✓ All modules imported successfully')
" || {
        error "Failed to import GenomicsToolkit modules!"
        exit 1
    }
    
    success "Installation validated"
}

# Function to handle signals
cleanup() {
    log "Received signal, cleaning up..."
    # Add any cleanup tasks here
    exit 0
}

# Function to display help
show_help() {
    cat << EOF
GenomicsToolkit Docker Container

Usage: docker run genomics-toolkit [OPTIONS] [COMMAND]

Commands:
  analyze-sequence    Analyze DNA/RNA sequences
  call-variants      Call variants from alignment data
  visualize          Generate visualizations
  pipeline           Run complete analysis pipeline
  config             Manage configuration
  validate           Validate files and system
  
Options:
  --help             Show this help message
  --version          Show version information
  --log-level LEVEL  Set logging level (DEBUG, INFO, WARNING, ERROR)
  
Environment Variables:
  GENOMICS_TOOLKIT_LOG_LEVEL     Logging level (default: INFO)
  GENOMICS_TOOLKIT_MAX_WORKERS   Number of worker processes (default: 4)
  GENOMICS_TOOLKIT_CONFIG        Configuration directory
  GENOMICS_TOOLKIT_DATA          Input data directory
  GENOMICS_TOOLKIT_OUTPUT        Output directory

Examples:
  # Analyze a FASTA file
  docker run -v /path/to/data:/app/data -v /path/to/output:/app/output \\
    genomics-toolkit analyze-sequence /app/data/sequences.fasta

  # Run complete pipeline
  docker run -v /path/to/data:/app/data -v /path/to/output:/app/output \\
    genomics-toolkit pipeline /app/data/input.fastq --reference /app/data/ref.fasta

  # Interactive shell
  docker run -it -v /path/to/data:/app/data genomics-toolkit bash

For more information, visit: https://github.com/yourusername/genomics-toolkit
EOF
}

# Trap signals
trap cleanup SIGTERM SIGINT

# Main execution
main() {
    log "Starting GenomicsToolkit container..."
    
    # Handle special cases
    case "${1:-}" in
        --help|-h)
            show_help
            exit 0
            ;;
        --version|-v)
            genomics-toolkit --version
            exit 0
            ;;
        bash|sh|/bin/bash|/bin/sh)
            log "Starting interactive shell..."
            exec "$@"
            ;;
    esac
    
    # Initialize container
    init_config
    setup_directories
    check_resources
    validate_installation
    
    success "GenomicsToolkit container ready!"
    
    # Execute the main command
    if [ $# -eq 0 ]; then
        # No command provided, show help
        genomics-toolkit --help
    else
        log "Executing: $*"
        exec "$@"
    fi
}

# Run main function with all arguments
main "$@"