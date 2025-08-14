"""
Sequence Analysis Module

This module provides comprehensive DNA/RNA sequence analysis capabilities including:
- FASTA/FASTQ parsing and validation
- Basic sequence statistics
- ORF detection
- Translation and reverse complement operations
- Quality score analysis
"""

import re
import logging
from typing import Dict, List, Tuple, Optional, Generator, Union
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
try:
    from Bio.SeqUtils import GC, molecular_weight
except ImportError:
    # Fallback if GC import fails
    def GC(seq):
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        return (gc_count / len(seq)) * 100 if seq else 0
    
    def molecular_weight(seq):
        return len(seq) * 650  # Approximate molecular weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class SequenceAnalyzer:
    """Main class for sequence analysis operations."""
    
    def __init__(self, log_level: str = "INFO"):
        """Initialize the SequenceAnalyzer."""
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        # Genetic code for translation
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        self.start_codons = {'ATG', 'GTG', 'TTG'}
        self.stop_codons = {'TAA', 'TAG', 'TGA'}

    def parse_fasta(self, file_path: Union[str, Path]) -> Generator[Tuple[str, str], None, None]:
        """
        Parse FASTA file and yield (header, sequence) tuples.
        
        Args:
            file_path: Path to FASTA file
            
        Yields:
            Tuple of (header, sequence)
        """
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                yield record.id, str(record.seq)
        except Exception as e:
            self.logger.error(f"Error parsing FASTA file {file_path}: {e}")
            raise

    def parse_fastq(self, file_path: Union[str, Path]) -> Generator[Tuple[str, str, str], None, None]:
        """
        Parse FASTQ file and yield (header, sequence, quality) tuples.
        
        Args:
            file_path: Path to FASTQ file
            
        Yields:
            Tuple of (header, sequence, quality)
        """
        try:
            for record in SeqIO.parse(file_path, "fastq"):
                yield record.id, str(record.seq), record.letter_annotations["phred_quality"]
        except Exception as e:
            self.logger.error(f"Error parsing FASTQ file {file_path}: {e}")
            raise

    def validate_sequence(self, sequence: str, seq_type: str = "DNA") -> bool:
        """
        Validate DNA/RNA sequence.
        
        Args:
            sequence: Sequence string
            seq_type: Type of sequence ("DNA", "RNA", or "protein")
            
        Returns:
            True if valid, False otherwise
        """
        sequence = sequence.upper().replace("U", "T") if seq_type == "RNA" else sequence.upper()
        
        if seq_type.upper() == "DNA":
            valid_chars = set("ATCGN-")
        elif seq_type.upper() == "RNA":
            valid_chars = set("AUCGN-")
        elif seq_type.upper() == "PROTEIN":
            valid_chars = set("ACDEFGHIKLMNPQRSTVWY*-")
        else:
            raise ValueError("seq_type must be 'DNA', 'RNA', or 'protein'")
            
        return all(char in valid_chars for char in sequence)

    def basic_stats(self, sequence: str) -> Dict[str, Union[int, float]]:
        """
        Calculate basic sequence statistics.
        
        Args:
            sequence: DNA/RNA sequence
            
        Returns:
            Dictionary with sequence statistics
        """
        sequence = sequence.upper()
        length = len(sequence)
        
        if length == 0:
            return {"length": 0, "gc_content": 0.0, "at_content": 0.0}
        
        composition = Counter(sequence)
        gc_count = composition.get('G', 0) + composition.get('C', 0)
        at_count = composition.get('A', 0) + composition.get('T', 0) + composition.get('U', 0)
        
        stats = {
            "length": length,
            "gc_content": (gc_count / length) * 100,
            "at_content": (at_count / length) * 100,
            "n_count": composition.get('N', 0),
            "composition": dict(composition)
        }
        
        return stats

    def reverse_complement(self, sequence: str) -> str:
        """
        Get reverse complement of DNA sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Reverse complement sequence
        """
        complement_map = str.maketrans("ATCGN", "TAGCN")
        return sequence.upper().translate(complement_map)[::-1]

    def transcribe(self, dna_sequence: str) -> str:
        """
        Transcribe DNA to RNA.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            RNA sequence
        """
        return dna_sequence.upper().replace("T", "U")

    def translate(self, sequence: str, frame: int = 0) -> str:
        """
        Translate DNA/RNA sequence to protein.
        
        Args:
            sequence: DNA/RNA sequence
            frame: Reading frame (0, 1, or 2)
            
        Returns:
            Protein sequence
        """
        sequence = sequence.upper().replace("U", "T")
        
        if frame < 0 or frame > 2:
            raise ValueError("Frame must be 0, 1, or 2")
            
        # Start from the specified frame
        seq_frame = sequence[frame:]
        
        # Translate in codons
        protein = ""
        for i in range(0, len(seq_frame) - 2, 3):
            codon = seq_frame[i:i+3]
            if len(codon) == 3:
                protein += self.genetic_code.get(codon, "X")
            
        return protein

    def find_orfs(self, sequence: str, min_length: int = 100) -> List[Dict]:
        """
        Find Open Reading Frames (ORFs) in all six frames.
        
        Args:
            sequence: DNA sequence
            min_length: Minimum ORF length in nucleotides
            
        Returns:
            List of ORF dictionaries
        """
        sequence = sequence.upper()
        orfs = []
        
        # Check all six reading frames
        for strand in [1, -1]:
            seq = sequence if strand == 1 else self.reverse_complement(sequence)
            
            for frame in range(3):
                frame_seq = seq[frame:]
                
                # Find ORFs in this frame
                i = 0
                while i < len(frame_seq) - 2:
                    codon = frame_seq[i:i+3]
                    
                    if len(codon) == 3 and codon in self.start_codons:
                        # Found start codon, look for stop codon
                        j = i + 3
                        while j < len(frame_seq) - 2:
                            stop_codon = frame_seq[j:j+3]
                            if len(stop_codon) == 3 and stop_codon in self.stop_codons:
                                # Found ORF
                                orf_length = j + 3 - i
                                if orf_length >= min_length:
                                    start_pos = i + frame if strand == 1 else len(sequence) - (i + frame + orf_length)
                                    
                                    orf = {
                                        "start": start_pos,
                                        "end": start_pos + orf_length,
                                        "length": orf_length,
                                        "strand": strand,
                                        "frame": frame + 1,
                                        "sequence": frame_seq[i:j+3],
                                        "protein": self.translate(frame_seq[i:j+3])
                                    }
                                    orfs.append(orf)
                                
                                i = j + 3
                                break
                            j += 3
                        else:
                            # No stop codon found
                            break
                    else:
                        i += 3
        
        return sorted(orfs, key=lambda x: x["length"], reverse=True)

    def analyze_quality_scores(self, quality_scores: List[int]) -> Dict[str, float]:
        """
        Analyze FASTQ quality scores.
        
        Args:
            quality_scores: List of Phred quality scores
            
        Returns:
            Quality statistics dictionary
        """
        if not quality_scores:
            return {}
            
        q_array = np.array(quality_scores)
        
        stats = {
            "mean_quality": float(np.mean(q_array)),
            "median_quality": float(np.median(q_array)),
            "min_quality": float(np.min(q_array)),
            "max_quality": float(np.max(q_array)),
            "q25": float(np.percentile(q_array, 25)),
            "q75": float(np.percentile(q_array, 75)),
            "bases_q20": int(np.sum(q_array >= 20)),
            "bases_q30": int(np.sum(q_array >= 30)),
            "percent_q20": float(np.sum(q_array >= 20) / len(q_array) * 100),
            "percent_q30": float(np.sum(q_array >= 30) / len(q_array) * 100)
        }
        
        return stats

    def calculate_n50(self, sequences: List[str]) -> int:
        """
        Calculate N50 statistic for a set of sequences.
        
        Args:
            sequences: List of sequence strings
            
        Returns:
            N50 value
        """
        if not sequences:
            return 0
            
        lengths = sorted([len(seq) for seq in sequences], reverse=True)
        total_length = sum(lengths)
        cumulative_length = 0
        
        for length in lengths:
            cumulative_length += length
            if cumulative_length >= total_length / 2:
                return length
                
        return 0

    def composition_analysis(self, sequence: str, window_size: int = 100) -> Dict[str, List[float]]:
        """
        Analyze sequence composition in sliding windows.
        
        Args:
            sequence: DNA/RNA sequence
            window_size: Size of sliding window
            
        Returns:
            Dictionary with composition data for each window
        """
        sequence = sequence.upper()
        windows = []
        
        for i in range(0, len(sequence) - window_size + 1, window_size):
            window = sequence[i:i + window_size]
            windows.append(window)
        
        composition_data = {
            "positions": [],
            "gc_content": [],
            "at_content": [],
            "n_content": []
        }
        
        for i, window in enumerate(windows):
            stats = self.basic_stats(window)
            composition_data["positions"].append(i * window_size)
            composition_data["gc_content"].append(stats["gc_content"])
            composition_data["at_content"].append(stats["at_content"])
            composition_data["n_content"].append((stats["n_count"] / len(window)) * 100)
        
        return composition_data