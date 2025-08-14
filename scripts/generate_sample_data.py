#!/usr/bin/env python3
"""
Generate sample datasets for GenomicsToolkit tutorials and testing.
"""

import random
import os
from pathlib import Path

def generate_sample_data():
    """Generate all sample datasets."""
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    
    # Generate sample FASTA
    generate_fasta(data_dir / "sample_sequences.fasta", num_seqs=5, seq_length=1000)
    
    # Generate sample FASTQ
    generate_fastq(data_dir / "sample_reads.fastq", num_reads=1000, read_length=150)
    
    # Generate reference genome
    generate_fasta(data_dir / "reference_genome.fasta", num_seqs=1, seq_length=5000)
    
    # Generate sample VCF
    generate_vcf(data_dir / "sample_variants.vcf", num_variants=50)
    
    print("Sample datasets generated successfully!")
    print(f"Files created in {data_dir}:")
    for file in data_dir.glob("*"):
        print(f"  - {file.name}")

def generate_fasta(filepath, num_seqs=5, seq_length=1000):
    """Generate sample FASTA file."""
    bases = ['A', 'T', 'C', 'G']
    
    with open(filepath, 'w') as f:
        for i in range(num_seqs):
            f.write(f">sequence_{i+1} Sample sequence {i+1}\n")
            
            # Generate sequence with some realistic patterns
            sequence = ""
            for j in range(seq_length):
                if j % 100 == 0 and j > 0:  # Add some ORFs
                    sequence += "ATG"  # Start codon
                    for _ in range(30):  # 30 codons
                        sequence += random.choice(['AAA', 'TTT', 'CCC', 'GGG', 'ATC', 'GCA'])
                    sequence += random.choice(['TAA', 'TAG', 'TGA'])  # Stop codon
                    j += 93  # Skip ahead
                else:
                    sequence += random.choice(bases)
            
            # Add some Ns occasionally
            sequence = ''.join([base if random.random() > 0.01 else 'N' for base in sequence])
            
            # Write in 80-character lines
            for k in range(0, len(sequence), 80):
                f.write(sequence[k:k+80] + '\n')

def generate_fastq(filepath, num_reads=1000, read_length=150):
    """Generate sample FASTQ file."""
    bases = ['A', 'T', 'C', 'G']
    
    with open(filepath, 'w') as f:
        for i in range(num_reads):
            # Header
            f.write(f"@read_{i+1}\n")
            
            # Sequence
            sequence = ''.join(random.choice(bases) for _ in range(read_length))
            f.write(sequence + '\n')
            
            # Plus line
            f.write('+\n')
            
            # Quality scores (simulate decreasing quality towards end)
            quality = ""
            for j in range(read_length):
                base_quality = max(20, 40 - (j // 10))  # Quality decreases
                quality += chr(33 + base_quality)
            f.write(quality + '\n')

def generate_vcf(filepath, num_variants=50):
    """Generate sample VCF file."""
    with open(filepath, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##reference=sample_reference\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Generate variants
        bases = ['A', 'T', 'C', 'G']
        chroms = ['chr1', 'chr2', 'chr3']
        
        positions = sorted(random.sample(range(1, 5000), num_variants))
        
        for i, pos in enumerate(positions):
            chrom = random.choice(chroms)
            ref = random.choice(bases)
            
            # 80% SNPs, 20% indels
            if random.random() < 0.8:
                alt = random.choice([b for b in bases if b != ref])
            else:
                if random.random() < 0.5:  # Insertion
                    alt = ref + random.choice(bases)
                else:  # Deletion
                    alt = ref
                    ref = ref + random.choice(bases)
            
            qual = random.uniform(20, 60)
            depth = random.randint(10, 100)
            af = random.uniform(0.1, 0.9)
            
            f.write(f"{chrom}\t{pos}\tvar_{i+1}\t{ref}\t{alt}\t{qual:.1f}\tPASS\t")
            f.write(f"DP={depth};AF={af:.3f}\n")

if __name__ == "__main__":
    generate_sample_data()