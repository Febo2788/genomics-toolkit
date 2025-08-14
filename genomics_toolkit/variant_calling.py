"""
Variant Calling Pipeline Module

This module provides comprehensive variant calling capabilities including:
- SNP and indel detection
- Variant quality filtering and annotation
- VCF file generation and parsing
- Population genetics statistics
- Hardy-Weinberg equilibrium testing
"""

import re
import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
from collections import defaultdict, Counter
from pathlib import Path
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    print("Warning: pysam not available. BAM file support disabled.")
from scipy.stats import chi2


class Variant:
    """Class representing a genetic variant."""
    
    def __init__(self, chrom: str, pos: int, ref: str, alt: str, quality: float = 0.0):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.info = {}
        self.genotypes = []
        self.annotations = {}
    
    def __str__(self):
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"
    
    def is_snp(self) -> bool:
        """Check if variant is a SNP."""
        return len(self.ref) == 1 and len(self.alt) == 1
    
    def is_indel(self) -> bool:
        """Check if variant is an indel."""
        return len(self.ref) != len(self.alt)
    
    def is_insertion(self) -> bool:
        """Check if variant is an insertion."""
        return len(self.alt) > len(self.ref)
    
    def is_deletion(self) -> bool:
        """Check if variant is a deletion."""
        return len(self.ref) > len(self.alt)


class VariantCaller:
    """Main class for variant calling and analysis operations."""
    
    def __init__(self, log_level: str = "INFO"):
        """Initialize the VariantCaller."""
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))
        self.variants = []
    
    def call_variants_from_alignment(self, bam_file: str, reference_fasta: str, 
                                   min_coverage: int = 10, min_quality: float = 20.0,
                                   min_variant_freq: float = 0.2) -> List[Variant]:
        """
        Call variants from BAM alignment file.
        
        Args:
            bam_file: Path to BAM file
            reference_fasta: Path to reference FASTA
            min_coverage: Minimum coverage depth
            min_quality: Minimum base quality
            min_variant_freq: Minimum variant frequency
            
        Returns:
            List of Variant objects
        """
        if not PYSAM_AVAILABLE:
            raise ImportError("pysam is required for BAM file processing. Please install pysam or use simulation mode.")
            
        variants = []
        
        try:
            # Open BAM file
            bamfile = pysam.AlignmentFile(bam_file, "rb")
            
            # Open reference FASTA
            ref_fasta = pysam.FastaFile(reference_fasta)
            
            # Iterate through each reference sequence
            for chrom in bamfile.references:
                self.logger.info(f"Processing chromosome {chrom}")
                
                # Get chromosome length
                chrom_length = bamfile.get_reference_length(chrom)
                
                # Pileup over the chromosome
                for pileupcolumn in bamfile.pileup(chrom, 0, chrom_length, 
                                                 min_base_quality=min_quality):
                    
                    if pileupcolumn.nsegments < min_coverage:
                        continue
                    
                    pos = pileupcolumn.reference_pos + 1  # 1-based position
                    ref_base = ref_fasta.fetch(chrom, pos-1, pos).upper()
                    
                    # Count bases at this position
                    base_counts = defaultdict(int)
                    total_reads = 0
                    
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                            base_counts[base] += 1
                            total_reads += 1
                    
                    # Check for variants
                    for alt_base, count in base_counts.items():
                        if alt_base != ref_base:
                            freq = count / total_reads
                            if freq >= min_variant_freq and count >= min_coverage:
                                # Calculate quality score (simplified)
                                quality = -10 * np.log10(max(1 - freq, 0.001))
                                
                                variant = Variant(chrom, pos, ref_base, alt_base, quality)
                                variant.info = {
                                    "DP": total_reads,
                                    "AD": count,
                                    "AF": freq
                                }
                                variants.append(variant)
            
            bamfile.close()
            ref_fasta.close()
            
        except Exception as e:
            self.logger.error(f"Error calling variants: {e}")
            
        self.variants = variants
        return variants
    
    def simulate_variants(self, reference_seq: str, num_snps: int = 100, 
                         num_indels: int = 20) -> List[Variant]:
        """
        Simulate variants for testing purposes.
        
        Args:
            reference_seq: Reference sequence
            num_snps: Number of SNPs to simulate
            num_indels: Number of indels to simulate
            
        Returns:
            List of simulated variants
        """
        variants = []
        seq_length = len(reference_seq)
        bases = ['A', 'T', 'C', 'G']
        
        # Generate SNPs
        snp_positions = np.random.choice(seq_length, size=num_snps, replace=False)
        for pos in snp_positions:
            ref_base = reference_seq[pos].upper()
            alt_base = np.random.choice([b for b in bases if b != ref_base])
            quality = np.random.uniform(20, 60)
            
            variant = Variant("chr1", pos + 1, ref_base, alt_base, quality)
            variant.info = {
                "DP": np.random.randint(20, 100),
                "AF": np.random.uniform(0.2, 0.8)
            }
            variants.append(variant)
        
        # Generate indels
        indel_positions = np.random.choice(
            [p for p in range(seq_length) if p not in snp_positions], 
            size=num_indels, replace=False
        )
        
        for pos in indel_positions:
            ref_base = reference_seq[pos].upper()
            
            if np.random.random() < 0.5:  # Insertion
                insert_length = np.random.randint(1, 5)
                insert_seq = ''.join(np.random.choice(bases, size=insert_length))
                alt_seq = ref_base + insert_seq
            else:  # Deletion
                del_length = np.random.randint(1, 5)
                end_pos = min(pos + del_length + 1, seq_length)
                ref_seq = reference_seq[pos:end_pos].upper()
                alt_seq = ref_base
                ref_base = ref_seq
            
            quality = np.random.uniform(20, 60)
            variant = Variant("chr1", pos + 1, ref_base, alt_seq, quality)
            variant.info = {
                "DP": np.random.randint(20, 100),
                "AF": np.random.uniform(0.2, 0.8)
            }
            variants.append(variant)
        
        self.variants = sorted(variants, key=lambda v: v.pos)
        return self.variants
    
    def filter_variants(self, min_quality: float = 30.0, min_depth: int = 10,
                       max_depth: int = 1000, min_allele_freq: float = 0.05) -> List[Variant]:
        """
        Filter variants based on quality criteria.
        
        Args:
            min_quality: Minimum variant quality
            min_depth: Minimum read depth
            max_depth: Maximum read depth
            min_allele_freq: Minimum allele frequency
            
        Returns:
            List of filtered variants
        """
        filtered_variants = []
        
        for variant in self.variants:
            # Quality filter
            if variant.quality < min_quality:
                continue
            
            # Depth filter
            depth = variant.info.get("DP", 0)
            if depth < min_depth or depth > max_depth:
                continue
            
            # Allele frequency filter
            af = variant.info.get("AF", 0)
            if af < min_allele_freq:
                continue
            
            filtered_variants.append(variant)
        
        return filtered_variants
    
    def annotate_variants(self, variants: List[Variant], annotation_db: Dict = None) -> List[Variant]:
        """
        Annotate variants with functional information.
        
        Args:
            variants: List of variants to annotate
            annotation_db: Database of annotations (simplified)
            
        Returns:
            List of annotated variants
        """
        for variant in variants:
            # Basic annotation based on variant type
            if variant.is_snp():
                variant.annotations["type"] = "SNP"
                variant.annotations["effect"] = self._predict_snp_effect(variant)
            elif variant.is_indel():
                variant.annotations["type"] = "INDEL"
                if variant.is_insertion():
                    variant.annotations["subtype"] = "insertion"
                    variant.annotations["effect"] = "frameshift" if (len(variant.alt) - len(variant.ref)) % 3 != 0 else "inframe_insertion"
                else:
                    variant.annotations["subtype"] = "deletion"
                    variant.annotations["effect"] = "frameshift" if (len(variant.ref) - len(variant.alt)) % 3 != 0 else "inframe_deletion"
            
            # Add population frequency (simulated)
            variant.annotations["population_freq"] = np.random.uniform(0.001, 0.5)
            
            # Clinical significance (simulated)
            effects = ["benign", "likely_benign", "uncertain", "likely_pathogenic", "pathogenic"]
            variant.annotations["clinical_significance"] = np.random.choice(effects, p=[0.4, 0.3, 0.2, 0.07, 0.03])
        
        return variants
    
    def _predict_snp_effect(self, variant: Variant) -> str:
        """
        Predict the effect of a SNP (simplified).
        
        Args:
            variant: SNP variant
            
        Returns:
            Predicted effect
        """
        # Simplified effect prediction
        effects = ["synonymous", "missense", "nonsense", "splice_site"]
        probabilities = [0.25, 0.65, 0.05, 0.05]
        return np.random.choice(effects, p=probabilities)
    
    def calculate_population_stats(self, variants: List[Variant]) -> Dict[str, float]:
        """
        Calculate population genetics statistics.
        
        Args:
            variants: List of variants with genotype information
            
        Returns:
            Dictionary of population statistics
        """
        stats = {}
        
        # Calculate allele frequencies
        all_allele_freqs = [v.info.get("AF", 0) for v in variants]
        stats["mean_allele_freq"] = np.mean(all_allele_freqs)
        stats["median_allele_freq"] = np.median(all_allele_freqs)
        
        # Count variant types
        snp_count = sum(1 for v in variants if v.is_snp())
        indel_count = sum(1 for v in variants if v.is_indel())
        
        stats["total_variants"] = len(variants)
        stats["snp_count"] = snp_count
        stats["indel_count"] = indel_count
        stats["snp_ratio"] = snp_count / len(variants) if variants else 0
        
        # Transition/transversion ratio for SNPs
        transitions = 0
        transversions = 0
        
        for variant in variants:
            if variant.is_snp():
                ref, alt = variant.ref, variant.alt
                if (ref in "AG" and alt in "AG") or (ref in "CT" and alt in "CT"):
                    transitions += 1
                else:
                    transversions += 1
        
        stats["transitions"] = transitions
        stats["transversions"] = transversions
        stats["ti_tv_ratio"] = transitions / transversions if transversions > 0 else 0
        
        return stats
    
    def hardy_weinberg_test(self, observed_genotypes: List[Tuple[str, str]]) -> Dict[str, float]:
        """
        Perform Hardy-Weinberg equilibrium test.
        
        Args:
            observed_genotypes: List of (allele1, allele2) tuples
            
        Returns:
            Dictionary with HWE test results
        """
        # Count genotypes
        genotype_counts = defaultdict(int)
        allele_counts = defaultdict(int)
        
        for allele1, allele2 in observed_genotypes:
            # Sort alleles for consistent genotype representation
            genotype = tuple(sorted([allele1, allele2]))
            genotype_counts[genotype] += 1
            
            allele_counts[allele1] += 1
            allele_counts[allele2] += 1
        
        total_individuals = len(observed_genotypes)
        total_alleles = total_individuals * 2
        
        # Calculate allele frequencies
        allele_freqs = {allele: count / total_alleles for allele, count in allele_counts.items()}
        
        # Calculate expected genotype frequencies under HWE
        alleles = list(allele_freqs.keys())
        expected_genotypes = {}
        
        for i, allele1 in enumerate(alleles):
            for j, allele2 in enumerate(alleles[i:], i):
                genotype = tuple(sorted([allele1, allele2]))
                if i == j:  # Homozygote
                    expected_freq = allele_freqs[allele1] ** 2
                else:  # Heterozygote
                    expected_freq = 2 * allele_freqs[allele1] * allele_freqs[allele2]
                
                expected_genotypes[genotype] = expected_freq * total_individuals
        
        # Chi-square test
        chi_square = 0
        for genotype in expected_genotypes:
            observed = genotype_counts.get(genotype, 0)
            expected = expected_genotypes[genotype]
            
            if expected > 0:
                chi_square += (observed - expected) ** 2 / expected
        
        # Degrees of freedom = number of genotypes - number of alleles
        df = len(expected_genotypes) - len(alleles)
        p_value = 1 - chi2.cdf(chi_square, df) if df > 0 else 1.0
        
        return {
            "chi_square": chi_square,
            "degrees_freedom": df,
            "p_value": p_value,
            "allele_frequencies": allele_freqs,
            "observed_genotypes": dict(genotype_counts),
            "expected_genotypes": expected_genotypes
        }
    
    def write_vcf(self, variants: List[Variant], output_file: str, reference_name: str = "Unknown"):
        """
        Write variants to VCF format file.
        
        Args:
            variants: List of variants to write
            output_file: Output VCF file path
            reference_name: Reference genome name
        """
        with open(output_file, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##reference={reference_name}\n")
            f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
            f.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
            f.write("##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Write variants
            for i, variant in enumerate(variants):
                chrom = variant.chrom
                pos = variant.pos
                variant_id = f"var_{i+1}"
                ref = variant.ref
                alt = variant.alt
                qual = f"{variant.quality:.1f}"
                filter_field = "PASS"
                
                # Build INFO field
                info_parts = []
                for key, value in variant.info.items():
                    info_parts.append(f"{key}={value}")
                info = ";".join(info_parts)
                
                f.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\n")
    
    def parse_vcf(self, vcf_file: str) -> List[Variant]:
        """
        Parse variants from VCF file.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            List of Variant objects
        """
        variants = []
        
        with open(vcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Parse variant line
                fields = line.split('\t')
                if len(fields) >= 8:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    qual = float(fields[5]) if fields[5] != '.' else 0.0
                    
                    variant = Variant(chrom, pos, ref, alt, qual)
                    
                    # Parse INFO field
                    info_field = fields[7]
                    if info_field != '.':
                        for info_part in info_field.split(';'):
                            if '=' in info_part:
                                key, value = info_part.split('=', 1)
                                try:
                                    # Try to convert to appropriate type
                                    if '.' in value:
                                        variant.info[key] = float(value)
                                    else:
                                        variant.info[key] = int(value)
                                except ValueError:
                                    variant.info[key] = value
                    
                    variants.append(variant)
        
        return variants