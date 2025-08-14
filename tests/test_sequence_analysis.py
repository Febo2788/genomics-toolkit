"""
Unit tests for sequence analysis module.
"""

import pytest
import tempfile
import os
from pathlib import Path

from genomics_toolkit.sequence_analysis import SequenceAnalyzer


class TestSequenceAnalyzer:
    """Test cases for SequenceAnalyzer class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = SequenceAnalyzer()
        self.test_dna = "ATCGATCGATCGATCG"
        self.test_dna_with_orfs = "ATGAAAAAATTTAAAGGGCCCTATGGGGATCTAG"
        self.test_rna = "AUCGAUCGAUCGAUCG"
        
    def test_validate_dna_sequence(self):
        """Test DNA sequence validation."""
        assert self.analyzer.validate_sequence("ATCG", "DNA") is True
        assert self.analyzer.validate_sequence("ATCGN", "DNA") is True
        assert self.analyzer.validate_sequence("ATCGX", "DNA") is False
        assert self.analyzer.validate_sequence("", "DNA") is True
        
    def test_validate_rna_sequence(self):
        """Test RNA sequence validation."""
        assert self.analyzer.validate_sequence("AUCG", "RNA") is True
        assert self.analyzer.validate_sequence("AUCGN", "RNA") is True
        assert self.analyzer.validate_sequence("AUCGX", "RNA") is False
        
    def test_validate_protein_sequence(self):
        """Test protein sequence validation."""
        assert self.analyzer.validate_sequence("ACDEFGHIKLMNPQRSTVWY", "protein") is True
        assert self.analyzer.validate_sequence("ACDEFGHIKLMNPQRSTVWY*", "protein") is True
        assert self.analyzer.validate_sequence("ACDEFGHIKLMNPQRSTVWYX", "protein") is False
        
    def test_basic_stats(self):
        """Test basic sequence statistics calculation."""
        stats = self.analyzer.basic_stats(self.test_dna)
        
        assert stats["length"] == 16
        assert stats["gc_content"] == 50.0  # 8 G/C out of 16
        assert stats["at_content"] == 50.0  # 8 A/T out of 16
        assert stats["n_count"] == 0
        assert "composition" in stats
        
    def test_basic_stats_empty_sequence(self):
        """Test basic stats with empty sequence."""
        stats = self.analyzer.basic_stats("")
        
        assert stats["length"] == 0
        assert stats["gc_content"] == 0.0
        assert stats["at_content"] == 0.0
        
    def test_reverse_complement(self):
        """Test reverse complement generation."""
        rev_comp = self.analyzer.reverse_complement("ATCG")
        assert rev_comp == "CGAT"
        
        rev_comp = self.analyzer.reverse_complement("AAATTTCCCGGG")
        assert rev_comp == "CCCGGGAAATTT"
        
    def test_transcribe(self):
        """Test DNA to RNA transcription."""
        rna = self.analyzer.transcribe("ATCG")
        assert rna == "AUCG"
        
        rna = self.analyzer.transcribe("TTTT")
        assert rna == "UUUU"
        
    def test_translate(self):
        """Test DNA/RNA to protein translation."""
        # Test with start codon
        protein = self.analyzer.translate("ATGAAATTT")  # Met-Lys-Phe
        assert protein == "MKF"
        
        # Test with different frame
        protein = self.analyzer.translate("AATGAAATTT", frame=1)  # Frame shift
        assert protein.startswith("MK")
        
        # Test with stop codon
        protein = self.analyzer.translate("ATGTAA")  # Met-Stop
        assert protein == "M*"
        
    def test_translate_invalid_frame(self):
        """Test translation with invalid frame."""
        with pytest.raises(ValueError):
            self.analyzer.translate("ATG", frame=3)
            
    def test_find_orfs(self):
        """Test ORF detection."""
        orfs = self.analyzer.find_orfs(self.test_dna_with_orfs, min_length=9)
        
        assert len(orfs) > 0
        
        # Check ORF properties
        for orf in orfs:
            assert "start" in orf
            assert "end" in orf
            assert "length" in orf
            assert "strand" in orf
            assert "frame" in orf
            assert "sequence" in orf
            assert "protein" in orf
            assert orf["length"] >= 9
            
    def test_find_orfs_no_orfs(self):
        """Test ORF detection with sequence containing no ORFs."""
        orfs = self.analyzer.find_orfs("AAAAAAAAAAAA", min_length=100)
        assert len(orfs) == 0
        
    def test_analyze_quality_scores(self):
        """Test quality score analysis."""
        quality_scores = [20, 25, 30, 35, 40, 15, 25, 30]
        stats = self.analyzer.analyze_quality_scores(quality_scores)
        
        assert "mean_quality" in stats
        assert "median_quality" in stats
        assert "min_quality" in stats
        assert "max_quality" in stats
        assert "percent_q20" in stats
        assert "percent_q30" in stats
        
        assert stats["min_quality"] == 15
        assert stats["max_quality"] == 40
        assert stats["bases_q20"] == 7  # 7 bases >= Q20
        assert stats["bases_q30"] == 4  # 4 bases >= Q30
        
    def test_analyze_quality_scores_empty(self):
        """Test quality score analysis with empty list."""
        stats = self.analyzer.analyze_quality_scores([])
        assert stats == {}
        
    def test_calculate_n50(self):
        """Test N50 calculation."""
        sequences = ["AAAA", "TTTTTT", "CCCCCCCC", "GG"]  # lengths: 4, 6, 8, 2
        n50 = self.analyzer.calculate_n50(sequences)
        
        # Total length: 20, half: 10
        # Sorted lengths: 8, 6, 4, 2
        # Cumulative: 8, 14 (>= 10), so N50 = 6
        assert n50 == 6
        
    def test_calculate_n50_empty(self):
        """Test N50 calculation with empty list."""
        n50 = self.analyzer.calculate_n50([])
        assert n50 == 0
        
    def test_composition_analysis(self):
        """Test composition analysis with sliding windows."""
        sequence = "A" * 50 + "G" * 50 + "C" * 50 + "T" * 50  # 200 bp
        composition = self.analyzer.composition_analysis(sequence, window_size=50)
        
        assert "positions" in composition
        assert "gc_content" in composition
        assert "at_content" in composition
        assert "n_content" in composition
        
        assert len(composition["positions"]) == 4  # 4 windows
        assert composition["gc_content"][0] == 0.0  # First window: all A
        assert composition["gc_content"][1] == 100.0  # Second window: all G
        
    def test_parse_fasta_file(self):
        """Test FASTA file parsing."""
        # Create temporary FASTA file
        fasta_content = ">seq1\nATCGATCG\n>seq2\nGCGCGCGC\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
            
        try:
            sequences = list(self.analyzer.parse_fasta(temp_path))
            
            assert len(sequences) == 2
            assert sequences[0] == ("seq1", "ATCGATCG")
            assert sequences[1] == ("seq2", "GCGCGCGC")
            
        finally:
            os.unlink(temp_path)
            
    def test_parse_fasta_invalid_file(self):
        """Test FASTA parsing with invalid file."""
        with pytest.raises(Exception):
            list(self.analyzer.parse_fasta("nonexistent_file.fasta"))
            
    def test_parse_fastq_file(self):
        """Test FASTQ file parsing."""
        # Create temporary FASTQ file
        fastq_content = "@seq1\nATCG\n+\nIIII\n@seq2\nGCGC\n+\nJJJJ\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(fastq_content)
            temp_path = f.name
            
        try:
            sequences = list(self.analyzer.parse_fastq(temp_path))
            
            assert len(sequences) == 2
            assert sequences[0][0] == "seq1"
            assert sequences[0][1] == "ATCG"
            assert len(sequences[0][2]) == 4  # Quality scores
            
        finally:
            os.unlink(temp_path)


class TestSequenceAnalyzerEdgeCases:
    """Test edge cases and error conditions."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = SequenceAnalyzer()
        
    def test_very_long_sequence(self):
        """Test with very long sequence."""
        long_sequence = "ATCG" * 10000  # 40,000 bp
        stats = self.analyzer.basic_stats(long_sequence)
        
        assert stats["length"] == 40000
        assert stats["gc_content"] == 50.0
        
    def test_sequence_with_ns(self):
        """Test sequence with N characters."""
        sequence_with_ns = "ATCGNNNNATCG"
        stats = self.analyzer.basic_stats(sequence_with_ns)
        
        assert stats["n_count"] == 4
        assert stats["length"] == 12
        
    def test_mixed_case_sequence(self):
        """Test sequence with mixed case."""
        mixed_case = "AtCgAtCg"
        stats = self.analyzer.basic_stats(mixed_case)
        
        # Should handle mixed case correctly
        assert stats["length"] == 8
        assert stats["gc_content"] == 50.0
        
    def test_orf_detection_edge_cases(self):
        """Test ORF detection edge cases."""
        # Sequence with start codon at the end
        sequence = "AAAAAATG"
        orfs = self.analyzer.find_orfs(sequence, min_length=3)
        assert len(orfs) == 0  # Not enough space for complete ORF
        
        # Sequence with stop codon but no start
        sequence = "AAAAAATAA"
        orfs = self.analyzer.find_orfs(sequence, min_length=3)
        assert len(orfs) == 0
        
    def test_translation_with_incomplete_codon(self):
        """Test translation with incomplete codon at end."""
        protein = self.analyzer.translate("ATGAAA")  # Complete codons
        assert protein == "MK"
        
        protein = self.analyzer.translate("ATGAAAA")  # Incomplete codon at end
        assert protein == "MK"  # Should ignore incomplete codon


if __name__ == "__main__":
    pytest.main([__file__])