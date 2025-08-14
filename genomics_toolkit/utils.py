"""
Utility Module

This module provides utility functions and classes for:
- File handling operations
- Quality control checks
- Database integrations
- Multi-threading support
- Configuration management
"""

import os
import gzip
import logging
import threading
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from typing import Dict, List, Tuple, Optional, Union, Any, Callable
from pathlib import Path
import requests
import time
import psutil
import yaml
import json
from functools import wraps
from tqdm import tqdm


class FileHandler:
    """Utility class for file operations."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def detect_file_format(self, file_path: Union[str, Path]) -> str:
        """
        Detect file format based on extension and content.
        
        Args:
            file_path: Path to the file
            
        Returns:
            File format string
        """
        file_path = Path(file_path)
        
        # Check extension
        if file_path.suffix.lower() in ['.fa', '.fasta']:
            return 'fasta'
        elif file_path.suffix.lower() in ['.fq', '.fastq']:
            return 'fastq'
        elif file_path.suffix.lower() == '.vcf':
            return 'vcf'
        elif file_path.suffix.lower() == '.bam':
            return 'bam'
        elif file_path.suffix.lower() == '.sam':
            return 'sam'
        elif file_path.suffix.lower() == '.bed':
            return 'bed'
        elif file_path.suffix.lower() == '.gff':
            return 'gff'
        
        # Try to detect from content
        try:
            with self._open_file(file_path) as f:
                first_line = f.readline().strip()
                
                if first_line.startswith('>'):
                    return 'fasta'
                elif first_line.startswith('@'):
                    return 'fastq'
                elif first_line.startswith('##fileformat=VCF'):
                    return 'vcf'
                elif first_line.startswith('@HD') or first_line.startswith('@SQ'):
                    return 'sam'
        except Exception as e:
            self.logger.warning(f"Could not detect format for {file_path}: {e}")
        
        return 'unknown'
    
    def _open_file(self, file_path: Union[str, Path], mode: str = 'rt'):
        """
        Open file with automatic gzip detection.
        
        Args:
            file_path: Path to file
            mode: File mode
            
        Returns:
            File handle
        """
        file_path = Path(file_path)
        
        if file_path.suffix.lower() == '.gz':
            return gzip.open(file_path, mode)
        else:
            return open(file_path, mode)
    
    def validate_file_integrity(self, file_path: Union[str, Path]) -> Dict[str, Any]:
        """
        Validate file integrity and provide basic statistics.
        
        Args:
            file_path: Path to file
            
        Returns:
            Dictionary with validation results
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            return {"valid": False, "error": "File does not exist"}
        
        try:
            file_size = file_path.stat().st_size
            file_format = self.detect_file_format(file_path)
            
            # Count lines for text files
            line_count = 0
            if file_format in ['fasta', 'fastq', 'vcf', 'sam', 'bed', 'gff']:
                with self._open_file(file_path) as f:
                    line_count = sum(1 for _ in f)
            
            return {
                "valid": True,
                "file_size": file_size,
                "file_format": file_format,
                "line_count": line_count,
                "readable": True
            }
            
        except Exception as e:
            return {"valid": False, "error": str(e)}
    
    def split_large_file(self, input_file: Union[str, Path], chunk_size: int = 10000,
                        output_dir: Optional[str] = None) -> List[Path]:
        """
        Split large file into smaller chunks.
        
        Args:
            input_file: Input file path
            chunk_size: Number of records per chunk
            output_dir: Output directory
            
        Returns:
            List of chunk file paths
        """
        input_path = Path(input_file)
        output_path = Path(output_dir) if output_dir else input_path.parent
        output_path.mkdir(exist_ok=True)
        
        file_format = self.detect_file_format(input_path)
        chunk_files = []
        
        if file_format == 'fasta':
            chunk_files = self._split_fasta(input_path, chunk_size, output_path)
        elif file_format == 'fastq':
            chunk_files = self._split_fastq(input_path, chunk_size, output_path)
        else:
            self.logger.warning(f"File splitting not implemented for {file_format}")
        
        return chunk_files
    
    def _split_fasta(self, input_file: Path, chunk_size: int, output_dir: Path) -> List[Path]:
        """Split FASTA file into chunks."""
        chunk_files = []
        chunk_num = 1
        record_count = 0
        
        current_chunk = None
        
        with self._open_file(input_file) as f:
            for line in f:
                if line.startswith('>'):
                    if record_count % chunk_size == 0:
                        if current_chunk:
                            current_chunk.close()
                        
                        chunk_file = output_dir / f"{input_file.stem}_chunk_{chunk_num}.fasta"
                        current_chunk = open(chunk_file, 'w')
                        chunk_files.append(chunk_file)
                        chunk_num += 1
                    
                    record_count += 1
                
                if current_chunk:
                    current_chunk.write(line)
        
        if current_chunk:
            current_chunk.close()
        
        return chunk_files
    
    def _split_fastq(self, input_file: Path, chunk_size: int, output_dir: Path) -> List[Path]:
        """Split FASTQ file into chunks."""
        chunk_files = []
        chunk_num = 1
        record_count = 0
        line_count = 0
        
        current_chunk = None
        
        with self._open_file(input_file) as f:
            for line in f:
                if line_count % 4 == 0:  # Header line
                    if record_count % chunk_size == 0:
                        if current_chunk:
                            current_chunk.close()
                        
                        chunk_file = output_dir / f"{input_file.stem}_chunk_{chunk_num}.fastq"
                        current_chunk = open(chunk_file, 'w')
                        chunk_files.append(chunk_file)
                        chunk_num += 1
                    
                    record_count += 1
                
                if current_chunk:
                    current_chunk.write(line)
                
                line_count += 1
        
        if current_chunk:
            current_chunk.close()
        
        return chunk_files


class QualityControl:
    """Quality control and validation utilities."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def check_sequence_quality(self, sequence: str, min_length: int = 50) -> Dict[str, Any]:
        """
        Perform quality checks on sequence.
        
        Args:
            sequence: DNA/RNA sequence
            min_length: Minimum acceptable length
            
        Returns:
            Quality check results
        """
        results = {
            "passed": True,
            "warnings": [],
            "errors": []
        }
        
        # Length check
        if len(sequence) < min_length:
            results["errors"].append(f"Sequence too short: {len(sequence)} < {min_length}")
            results["passed"] = False
        
        # N content check
        n_count = sequence.upper().count('N')
        n_percentage = (n_count / len(sequence)) * 100 if sequence else 0
        
        if n_percentage > 50:
            results["errors"].append(f"High N content: {n_percentage:.1f}%")
            results["passed"] = False
        elif n_percentage > 10:
            results["warnings"].append(f"Elevated N content: {n_percentage:.1f}%")
        
        # Character validation
        valid_chars = set('ATCGN')
        invalid_chars = set(sequence.upper()) - valid_chars
        if invalid_chars:
            results["warnings"].append(f"Invalid characters found: {', '.join(invalid_chars)}")
        
        return results
    
    def check_system_resources(self) -> Dict[str, Any]:
        """
        Check available system resources.
        
        Returns:
            System resource information
        """
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage('/')
        
        return {
            "cpu_count": psutil.cpu_count(),
            "memory_total": memory.total,
            "memory_available": memory.available,
            "memory_percent": memory.percent,
            "disk_total": disk.total,
            "disk_free": disk.free,
            "disk_percent": disk.percent
        }
    
    def estimate_processing_time(self, file_size: int, operation_type: str = "sequence_analysis") -> float:
        """
        Estimate processing time based on file size.
        
        Args:
            file_size: File size in bytes
            operation_type: Type of operation
            
        Returns:
            Estimated time in seconds
        """
        # Rough estimates based on operation type
        rates = {
            "sequence_analysis": 1e6,  # 1MB per second
            "variant_calling": 5e5,    # 500KB per second
            "alignment": 1e5           # 100KB per second
        }
        
        rate = rates.get(operation_type, 1e6)
        return file_size / rate


class DatabaseConnector:
    """Connector for public biological databases."""
    
    def __init__(self, cache_dir: Optional[str] = None):
        self.logger = logging.getLogger(__name__)
        self.cache_dir = Path(cache_dir) if cache_dir else Path.home() / '.genomics_toolkit_cache'
        self.cache_dir.mkdir(exist_ok=True)
        
        # Base URLs for common databases
        self.base_urls = {
            'ncbi': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
            'ensembl': 'https://rest.ensembl.org/',
            'uniprot': 'https://www.uniprot.org/uniprot/'
        }
    
    def fetch_sequence(self, accession: str, database: str = 'ncbi') -> Optional[str]:
        """
        Fetch sequence from public database.
        
        Args:
            accession: Sequence accession number
            database: Database name ('ncbi', 'ensembl', 'uniprot')
            
        Returns:
            Sequence string or None
        """
        cache_file = self.cache_dir / f"{database}_{accession}.fasta"
        
        # Check cache first
        if cache_file.exists():
            with open(cache_file, 'r') as f:
                content = f.read()
                # Extract sequence from FASTA
                lines = content.split('\n')
                sequence = ''.join(line for line in lines if not line.startswith('>'))
                return sequence
        
        try:
            if database == 'ncbi':
                sequence = self._fetch_ncbi_sequence(accession)
            elif database == 'ensembl':
                sequence = self._fetch_ensembl_sequence(accession)
            else:
                self.logger.error(f"Unsupported database: {database}")
                return None
            
            # Cache the result
            if sequence:
                with open(cache_file, 'w') as f:
                    f.write(f">{accession}\n{sequence}\n")
            
            return sequence
            
        except Exception as e:
            self.logger.error(f"Error fetching sequence {accession} from {database}: {e}")
            return None
    
    def _fetch_ncbi_sequence(self, accession: str) -> Optional[str]:
        """Fetch sequence from NCBI."""
        url = f"{self.base_urls['ncbi']}efetch.fcgi"
        params = {
            'db': 'nucleotide',
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text'
        }
        
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        
        # Parse FASTA
        lines = response.text.split('\n')
        sequence = ''.join(line for line in lines if not line.startswith('>'))
        return sequence.replace('\n', '').replace('\r', '')
    
    def _fetch_ensembl_sequence(self, gene_id: str) -> Optional[str]:
        """Fetch sequence from Ensembl."""
        url = f"{self.base_urls['ensembl']}sequence/id/{gene_id}"
        headers = {'Content-Type': 'application/json'}
        
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        return data.get('seq')


class ParallelProcessor:
    """Utility for parallel processing of large datasets."""
    
    def __init__(self, max_workers: Optional[int] = None):
        self.max_workers = max_workers or mp.cpu_count()
        self.logger = logging.getLogger(__name__)
    
    def process_sequences_parallel(self, sequences: List[str], 
                                 func: Callable, 
                                 use_processes: bool = True,
                                 chunk_size: int = 100) -> List[Any]:
        """
        Process sequences in parallel.
        
        Args:
            sequences: List of sequences to process
            func: Function to apply to each sequence
            use_processes: Use processes instead of threads
            chunk_size: Size of chunks for processing
            
        Returns:
            List of results
        """
        if use_processes:
            executor_class = ProcessPoolExecutor
        else:
            executor_class = ThreadPoolExecutor
        
        results = []
        
        with executor_class(max_workers=self.max_workers) as executor:
            # Process in chunks to avoid memory issues
            for i in tqdm(range(0, len(sequences), chunk_size), desc="Processing chunks"):
                chunk = sequences[i:i+chunk_size]
                futures = [executor.submit(func, seq) for seq in chunk]
                
                for future in futures:
                    try:
                        result = future.result(timeout=300)  # 5 minute timeout
                        results.append(result)
                    except Exception as e:
                        self.logger.error(f"Error processing sequence: {e}")
                        results.append(None)
        
        return results


def timing_decorator(func: Callable) -> Callable:
    """Decorator to measure function execution time."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        
        logger = logging.getLogger(func.__module__)
        logger.info(f"{func.__name__} executed in {end_time - start_time:.2f} seconds")
        
        return result
    return wrapper


def retry_decorator(max_retries: int = 3, delay: float = 1.0):
    """Decorator to retry function calls on failure."""
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None
            
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_exception = e
                    if attempt < max_retries - 1:
                        time.sleep(delay * (2 ** attempt))  # Exponential backoff
                    
            raise last_exception
        return wrapper
    return decorator


class ConfigManager:
    """Configuration management utility."""
    
    def __init__(self, config_file: Optional[str] = None):
        self.config_file = Path(config_file) if config_file else Path.home() / '.genomics_toolkit_config.yaml'
        self.config = self._load_config()
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file."""
        if self.config_file.exists():
            with open(self.config_file, 'r') as f:
                return yaml.safe_load(f) or {}
        else:
            return self._create_default_config()
    
    def _create_default_config(self) -> Dict[str, Any]:
        """Create default configuration."""
        default_config = {
            'logging': {
                'level': 'INFO',
                'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            },
            'processing': {
                'max_workers': mp.cpu_count(),
                'chunk_size': 1000,
                'memory_limit': '8GB'
            },
            'quality_control': {
                'min_sequence_length': 50,
                'max_n_content': 10,
                'min_quality_score': 20
            },
            'databases': {
                'cache_dir': str(Path.home() / '.genomics_toolkit_cache'),
                'ncbi_email': '',
                'request_delay': 0.5
            }
        }
        
        # Save default config
        with open(self.config_file, 'w') as f:
            yaml.dump(default_config, f, default_flow_style=False)
        
        return default_config
    
    def get(self, key: str, default=None):
        """Get configuration value."""
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def set(self, key: str, value: Any):
        """Set configuration value."""
        keys = key.split('.')
        config = self.config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
        
        # Save updated config
        with open(self.config_file, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False)