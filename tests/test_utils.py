"""
Unit tests for utility modules.
"""

import pytest
import tempfile
import os
import json
from pathlib import Path
from unittest.mock import patch, MagicMock

from genomics_toolkit.utils import (
    FileHandler, QualityControl, DatabaseConnector, 
    ParallelProcessor, ConfigManager, timing_decorator, retry_decorator
)


class TestFileHandler:
    """Test cases for FileHandler class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.handler = FileHandler()
        
    def test_detect_fasta_format(self):
        """Test FASTA format detection."""
        # Test by extension
        assert self.handler.detect_file_format("test.fasta") == "fasta"
        assert self.handler.detect_file_format("test.fa") == "fasta"
        
        # Test by content
        fasta_content = ">seq1\nATCG\n>seq2\nGCGC\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
            
        try:
            assert self.handler.detect_file_format(temp_path) == "fasta"
        finally:
            os.unlink(temp_path)
            
    def test_detect_fastq_format(self):
        """Test FASTQ format detection."""
        # Test by extension
        assert self.handler.detect_file_format("test.fastq") == "fastq"
        assert self.handler.detect_file_format("test.fq") == "fastq"
        
        # Test by content
        fastq_content = "@seq1\nATCG\n+\nIIII\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(fastq_content)
            temp_path = f.name
            
        try:
            assert self.handler.detect_file_format(temp_path) == "fastq"
        finally:
            os.unlink(temp_path)
            
    def test_detect_vcf_format(self):
        """Test VCF format detection."""
        assert self.handler.detect_file_format("test.vcf") == "vcf"
        
        # Test by content
        vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(vcf_content)
            temp_path = f.name
            
        try:
            assert self.handler.detect_file_format(temp_path) == "vcf"
        finally:
            os.unlink(temp_path)
            
    def test_detect_unknown_format(self):
        """Test unknown format detection."""
        assert self.handler.detect_file_format("test.xyz") == "unknown"
        
    def test_validate_file_integrity(self):
        """Test file integrity validation."""
        # Create a test file
        content = "ATCGATCG\nGCGCGCGC\n"
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(content)
            temp_path = f.name
            
        try:
            result = self.handler.validate_file_integrity(temp_path)
            
            assert result["valid"] is True
            assert result["file_size"] > 0
            assert result["line_count"] == 2
            assert result["readable"] is True
            
        finally:
            os.unlink(temp_path)
            
    def test_validate_nonexistent_file(self):
        """Test validation of nonexistent file."""
        result = self.handler.validate_file_integrity("nonexistent_file.txt")
        
        assert result["valid"] is False
        assert "does not exist" in result["error"]
        
    def test_split_fasta_file(self):
        """Test FASTA file splitting."""
        # Create a test FASTA file with multiple sequences
        fasta_content = ">seq1\nATCG\n>seq2\nGCGC\n>seq3\nTTTT\n>seq4\nAAAA\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
            
        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                chunk_files = self.handler.split_large_file(
                    temp_path, chunk_size=2, output_dir=temp_dir
                )
                
                assert len(chunk_files) == 2  # 4 sequences / 2 per chunk = 2 chunks
                
                # Check that chunk files exist and contain sequences
                for chunk_file in chunk_files:
                    assert chunk_file.exists()
                    with open(chunk_file, 'r') as f:
                        content = f.read()
                        assert content.count('>') == 2  # 2 sequences per chunk
                        
            finally:
                os.unlink(temp_path)


class TestQualityControl:
    """Test cases for QualityControl class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.qc = QualityControl()
        
    def test_check_sequence_quality_good(self):
        """Test quality check with good sequence."""
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # 52 bp
        result = self.qc.check_sequence_quality(sequence, min_length=50)
        
        assert result["passed"] is True
        assert len(result["warnings"]) == 0
        assert len(result["errors"]) == 0
        
    def test_check_sequence_quality_too_short(self):
        """Test quality check with short sequence."""
        sequence = "ATCG"  # 4 bp
        result = self.qc.check_sequence_quality(sequence, min_length=50)
        
        assert result["passed"] is False
        assert len(result["errors"]) > 0
        assert "too short" in result["errors"][0]
        
    def test_check_sequence_quality_high_n_content(self):
        """Test quality check with high N content."""
        sequence = "N" * 60 + "ATCG" * 10  # 60% N content
        result = self.qc.check_sequence_quality(sequence)
        
        assert result["passed"] is False
        assert any("High N content" in error for error in result["errors"])
        
    def test_check_sequence_quality_moderate_n_content(self):
        """Test quality check with moderate N content."""
        sequence = "N" * 15 + "ATCG" * 35  # 15% N content
        result = self.qc.check_sequence_quality(sequence)
        
        assert result["passed"] is True
        assert any("Elevated N content" in warning for warning in result["warnings"])
        
    def test_check_sequence_quality_invalid_characters(self):
        """Test quality check with invalid characters."""
        sequence = "ATCGXYZ"
        result = self.qc.check_sequence_quality(sequence)
        
        assert any("Invalid characters" in warning for warning in result["warnings"])
        
    def test_check_system_resources(self):
        """Test system resource checking."""
        resources = self.qc.check_system_resources()
        
        assert "cpu_count" in resources
        assert "memory_total" in resources
        assert "memory_available" in resources
        assert "memory_percent" in resources
        assert "disk_total" in resources
        assert "disk_free" in resources
        assert "disk_percent" in resources
        
        assert resources["cpu_count"] > 0
        assert resources["memory_total"] > 0
        assert resources["disk_total"] > 0
        
    def test_estimate_processing_time(self):
        """Test processing time estimation."""
        file_size = 1000000  # 1MB
        
        # Test different operation types
        time_seq = self.qc.estimate_processing_time(file_size, "sequence_analysis")
        time_var = self.qc.estimate_processing_time(file_size, "variant_calling")
        time_align = self.qc.estimate_processing_time(file_size, "alignment")
        
        assert time_seq > 0
        assert time_var > 0
        assert time_align > 0
        
        # Alignment should take longer than sequence analysis
        assert time_align > time_seq


class TestDatabaseConnector:
    """Test cases for DatabaseConnector class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        with tempfile.TemporaryDirectory() as temp_dir:
            self.connector = DatabaseConnector(cache_dir=temp_dir)
            
    @patch('genomics_toolkit.utils.requests.get')
    def test_fetch_ncbi_sequence(self, mock_get):
        """Test NCBI sequence fetching."""
        # Mock successful response
        mock_response = MagicMock()
        mock_response.text = ">test_accession\nATCGATCG\n"
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        sequence = self.connector.fetch_sequence("test_accession", "ncbi")
        
        assert sequence == "ATCGATCG"
        mock_get.assert_called_once()
        
    @patch('genomics_toolkit.utils.requests.get')
    def test_fetch_ensembl_sequence(self, mock_get):
        """Test Ensembl sequence fetching."""
        # Mock successful response
        mock_response = MagicMock()
        mock_response.json.return_value = {"seq": "ATCGATCG"}
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        sequence = self.connector.fetch_sequence("test_gene", "ensembl")
        
        assert sequence == "ATCGATCG"
        mock_get.assert_called_once()
        
    def test_fetch_sequence_unsupported_database(self):
        """Test fetching from unsupported database."""
        sequence = self.connector.fetch_sequence("test", "unsupported_db")
        assert sequence is None


class TestParallelProcessor:
    """Test cases for ParallelProcessor class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.processor = ParallelProcessor(max_workers=2)
        
    def test_process_sequences_parallel(self):
        """Test parallel sequence processing."""
        sequences = ["ATCG", "GCTA", "TTTT", "AAAA"]
        
        # Simple function to get sequence length
        def get_length(seq):
            return len(seq)
            
        results = self.processor.process_sequences_parallel(
            sequences, get_length, use_processes=False, chunk_size=2
        )
        
        assert len(results) == 4
        assert all(r == 4 for r in results)  # All sequences have length 4
        
    def test_process_sequences_with_exception(self):
        """Test parallel processing with exceptions."""
        sequences = ["ATCG", "GCTA"]
        
        # Function that raises exception
        def failing_function(seq):
            if seq == "GCTA":
                raise ValueError("Test exception")
            return len(seq)
            
        results = self.processor.process_sequences_parallel(
            sequences, failing_function, use_processes=False
        )
        
        assert len(results) == 2
        assert results[0] == 4  # First sequence processed successfully
        assert results[1] is None  # Second sequence failed


class TestConfigManager:
    """Test cases for ConfigManager class."""
    
    def test_config_manager_creation(self):
        """Test configuration manager creation."""
        with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as f:
            config_path = f.name
            
        try:
            config_manager = ConfigManager(config_path)
            
            # Should create default config
            assert config_manager.config is not None
            assert "logging" in config_manager.config
            assert "processing" in config_manager.config
            
        finally:
            if os.path.exists(config_path):
                os.unlink(config_path)
                
    def test_config_get_set(self):
        """Test getting and setting configuration values."""
        with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as f:
            config_path = f.name
            
        try:
            config_manager = ConfigManager(config_path)
            
            # Test setting and getting values
            config_manager.set("test.value", "test_data")
            assert config_manager.get("test.value") == "test_data"
            
            # Test getting non-existent value
            assert config_manager.get("non.existent", "default") == "default"
            
            # Test nested values
            assert config_manager.get("logging.level") == "INFO"
            
        finally:
            if os.path.exists(config_path):
                os.unlink(config_path)


class TestDecorators:
    """Test cases for decorator functions."""
    
    def test_timing_decorator(self):
        """Test timing decorator."""
        @timing_decorator
        def test_function():
            return "result"
            
        result = test_function()
        assert result == "result"
        
    def test_retry_decorator(self):
        """Test retry decorator."""
        call_count = 0
        
        @retry_decorator(max_retries=3, delay=0.01)
        def failing_function():
            nonlocal call_count
            call_count += 1
            if call_count < 3:
                raise ValueError("Test exception")
            return "success"
            
        result = failing_function()
        assert result == "success"
        assert call_count == 3
        
    def test_retry_decorator_max_retries_exceeded(self):
        """Test retry decorator when max retries exceeded."""
        @retry_decorator(max_retries=2, delay=0.01)
        def always_failing_function():
            raise ValueError("Always fails")
            
        with pytest.raises(ValueError):
            always_failing_function()


if __name__ == "__main__":
    pytest.main([__file__])