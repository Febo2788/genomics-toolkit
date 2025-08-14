"""
Unit tests for variant calling module.
"""

import pytest
import tempfile
import os
from pathlib import Path

from genomics_toolkit.variant_calling import VariantCaller, Variant


class TestVariant:
    """Test cases for Variant class."""
    
    def test_variant_creation(self):
        """Test variant object creation."""
        variant = Variant("chr1", 100, "A", "T", 30.0)
        
        assert variant.chrom == "chr1"
        assert variant.pos == 100
        assert variant.ref == "A"
        assert variant.alt == "T"
        assert variant.quality == 30.0
        
    def test_variant_string_representation(self):
        """Test variant string representation."""
        variant = Variant("chr1", 100, "A", "T")
        assert str(variant) == "chr1:100:A>T"
        
    def test_is_snp(self):
        """Test SNP detection."""
        snp = Variant("chr1", 100, "A", "T")
        assert snp.is_snp() is True
        
        insertion = Variant("chr1", 100, "A", "AT")
        assert insertion.is_snp() is False
        
        deletion = Variant("chr1", 100, "AT", "A")
        assert deletion.is_snp() is False
        
    def test_is_indel(self):
        """Test indel detection."""
        snp = Variant("chr1", 100, "A", "T")
        assert snp.is_indel() is False
        
        insertion = Variant("chr1", 100, "A", "AT")
        assert insertion.is_indel() is True
        
        deletion = Variant("chr1", 100, "AT", "A")
        assert deletion.is_indel() is True
        
    def test_is_insertion(self):
        """Test insertion detection."""
        insertion = Variant("chr1", 100, "A", "AT")
        assert insertion.is_insertion() is True
        
        deletion = Variant("chr1", 100, "AT", "A")
        assert insertion.is_insertion() is True
        assert deletion.is_insertion() is False
        
        snp = Variant("chr1", 100, "A", "T")
        assert snp.is_insertion() is False
        
    def test_is_deletion(self):
        """Test deletion detection."""
        deletion = Variant("chr1", 100, "AT", "A")
        assert deletion.is_deletion() is True
        
        insertion = Variant("chr1", 100, "A", "AT")
        assert insertion.is_deletion() is False
        
        snp = Variant("chr1", 100, "A", "T")
        assert snp.is_deletion() is False


class TestVariantCaller:
    """Test cases for VariantCaller class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.caller = VariantCaller()
        self.test_sequence = "ATCGATCGATCGATCGATCGATCG" * 100  # 2400 bp
        
    def test_variant_caller_creation(self):
        """Test variant caller initialization."""
        caller = VariantCaller()
        assert caller.variants == []
        
    def test_simulate_variants(self):
        """Test variant simulation."""
        variants = self.caller.simulate_variants(
            self.test_sequence, num_snps=10, num_indels=5
        )
        
        assert len(variants) == 15  # 10 SNPs + 5 indels
        
        # Check that we have both SNPs and indels
        snp_count = sum(1 for v in variants if v.is_snp())
        indel_count = sum(1 for v in variants if v.is_indel())
        
        assert snp_count == 10
        assert indel_count == 5
        
        # Check variant properties
        for variant in variants:
            assert variant.chrom == "chr1"
            assert 1 <= variant.pos <= len(self.test_sequence)
            assert variant.quality > 0
            assert "DP" in variant.info
            assert "AF" in variant.info
            
    def test_simulate_variants_empty_sequence(self):
        """Test variant simulation with empty sequence."""
        with pytest.raises(Exception):
            self.caller.simulate_variants("", num_snps=1, num_indels=1)
            
    def test_filter_variants(self):
        """Test variant filtering."""
        # Create test variants
        variants = [
            Variant("chr1", 100, "A", "T", 40.0),  # High quality
            Variant("chr1", 200, "C", "G", 15.0),  # Low quality
            Variant("chr1", 300, "G", "C", 35.0),  # Medium quality
        ]
        
        # Add info data
        variants[0].info = {"DP": 50, "AF": 0.3}
        variants[1].info = {"DP": 5, "AF": 0.1}   # Low depth
        variants[2].info = {"DP": 30, "AF": 0.02}  # Low allele freq
        
        self.caller.variants = variants
        
        # Filter with strict criteria
        filtered = self.caller.filter_variants(
            min_quality=30.0,
            min_depth=10,
            min_allele_freq=0.05
        )
        
        assert len(filtered) == 1  # Only first variant should pass
        assert filtered[0].quality == 40.0
        
    def test_annotate_variants(self):
        """Test variant annotation."""
        variants = [
            Variant("chr1", 100, "A", "T"),  # SNP
            Variant("chr1", 200, "C", "CT"), # Insertion
            Variant("chr1", 300, "AG", "A"), # Deletion
        ]
        
        annotated = self.caller.annotate_variants(variants)
        
        assert len(annotated) == 3
        
        # Check annotations
        for variant in annotated:
            assert "type" in variant.annotations
            assert "effect" in variant.annotations
            assert "population_freq" in variant.annotations
            assert "clinical_significance" in variant.annotations
            
        # Check specific annotations
        assert annotated[0].annotations["type"] == "SNP"
        assert annotated[1].annotations["type"] == "INDEL"
        assert annotated[1].annotations["subtype"] == "insertion"
        assert annotated[2].annotations["type"] == "INDEL"
        assert annotated[2].annotations["subtype"] == "deletion"
        
    def test_calculate_population_stats(self):
        """Test population statistics calculation."""
        variants = []
        
        # Create SNPs with different characteristics
        for i in range(10):
            variant = Variant("chr1", i*100, "A", "T")
            variant.info = {"AF": 0.1 + i*0.05}
            variants.append(variant)
            
        # Create some indels
        for i in range(5):
            variant = Variant("chr1", 1000 + i*100, "A", "AT")
            variant.info = {"AF": 0.2 + i*0.1}
            variants.append(variant)
            
        stats = self.caller.calculate_population_stats(variants)
        
        assert "total_variants" in stats
        assert "snp_count" in stats
        assert "indel_count" in stats
        assert "snp_ratio" in stats
        assert "mean_allele_freq" in stats
        assert "transitions" in stats
        assert "transversions" in stats
        assert "ti_tv_ratio" in stats
        
        assert stats["total_variants"] == 15
        assert stats["snp_count"] == 10
        assert stats["indel_count"] == 5
        assert stats["snp_ratio"] == 10/15
        
    def test_hardy_weinberg_test(self):
        """Test Hardy-Weinberg equilibrium test."""
        # Create genotype data
        genotypes = [
            ("A", "A"), ("A", "A"), ("A", "A"), ("A", "A"),  # 4 AA
            ("A", "T"), ("A", "T"), ("A", "T"), ("A", "T"),  # 4 AT
            ("T", "T"), ("T", "T")                            # 2 TT
        ]
        
        hwe_result = self.caller.hardy_weinberg_test(genotypes)
        
        assert "chi_square" in hwe_result
        assert "degrees_freedom" in hwe_result
        assert "p_value" in hwe_result
        assert "allele_frequencies" in hwe_result
        assert "observed_genotypes" in hwe_result
        assert "expected_genotypes" in hwe_result
        
        # Check allele frequencies
        allele_freqs = hwe_result["allele_frequencies"]
        assert "A" in allele_freqs
        assert "T" in allele_freqs
        assert abs(allele_freqs["A"] - 0.6) < 0.01  # 12/20 = 0.6
        assert abs(allele_freqs["T"] - 0.4) < 0.01  # 8/20 = 0.4
        
    def test_write_and_parse_vcf(self):
        """Test VCF file writing and parsing."""
        # Create test variants
        variants = [
            Variant("chr1", 100, "A", "T", 30.0),
            Variant("chr1", 200, "C", "G", 25.0),
            Variant("chr2", 300, "G", "C", 40.0),
        ]
        
        # Add info data
        for i, variant in enumerate(variants):
            variant.info = {
                "DP": 20 + i*10,
                "AF": 0.1 + i*0.1,
                "AD": 5 + i*2
            }
        
        # Write to temporary VCF file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            temp_vcf_path = f.name
            
        try:
            self.caller.write_vcf(variants, temp_vcf_path)
            
            # Check file was created
            assert os.path.exists(temp_vcf_path)
            
            # Parse the VCF file
            parsed_variants = self.caller.parse_vcf(temp_vcf_path)
            
            assert len(parsed_variants) == 3
            
            # Check parsed variants
            for i, variant in enumerate(parsed_variants):
                assert variant.chrom == variants[i].chrom
                assert variant.pos == variants[i].pos
                assert variant.ref == variants[i].ref
                assert variant.alt == variants[i].alt
                assert abs(variant.quality - variants[i].quality) < 0.1
                
        finally:
            if os.path.exists(temp_vcf_path):
                os.unlink(temp_vcf_path)
                
    def test_parse_vcf_invalid_file(self):
        """Test VCF parsing with invalid file."""
        with pytest.raises(Exception):
            self.caller.parse_vcf("nonexistent_file.vcf")
            
    def test_predict_snp_effect(self):
        """Test SNP effect prediction."""
        variant = Variant("chr1", 100, "A", "T")
        effect = self.caller._predict_snp_effect(variant)
        
        # Should return one of the predefined effects
        valid_effects = ["synonymous", "missense", "nonsense", "splice_site"]
        assert effect in valid_effects


class TestVariantCallerIntegration:
    """Integration tests for variant calling workflow."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.caller = VariantCaller()
        
    def test_complete_variant_calling_workflow(self):
        """Test complete variant calling workflow."""
        reference_sequence = "ATCGATCG" * 1000  # 8000 bp
        
        # Step 1: Simulate variants
        variants = self.caller.simulate_variants(
            reference_sequence, num_snps=50, num_indels=10
        )
        assert len(variants) == 60
        
        # Step 2: Filter variants
        filtered_variants = self.caller.filter_variants(
            min_quality=25.0,
            min_depth=15
        )
        assert len(filtered_variants) <= len(variants)
        
        # Step 3: Annotate variants
        annotated_variants = self.caller.annotate_variants(filtered_variants)
        assert len(annotated_variants) == len(filtered_variants)
        
        # Step 4: Calculate statistics
        stats = self.caller.calculate_population_stats(annotated_variants)
        assert stats["total_variants"] == len(annotated_variants)
        
        # Step 5: Write to VCF
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
            temp_path = f.name
            
        try:
            self.caller.write_vcf(annotated_variants, temp_path)
            
            # Step 6: Parse VCF back
            parsed_variants = self.caller.parse_vcf(temp_path)
            assert len(parsed_variants) == len(annotated_variants)
            
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)


if __name__ == "__main__":
    pytest.main([__file__])