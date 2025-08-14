"""
GenomicsToolkit: A comprehensive bioinformatics pipeline for DNA/RNA sequence analysis and variant calling.

This package provides tools for:
- FASTA/FASTQ file parsing and validation
- Sequence analysis and statistics
- Open Reading Frame (ORF) detection
- Variant calling and analysis
- Visualization and reporting
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

from .sequence_analysis import SequenceAnalyzer
from .variant_calling import VariantCaller
from .visualization import Visualizer
from .utils import FileHandler, QualityControl

try:
    from .gui import GenomicsToolkitGUI
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False

__all__ = [
    "SequenceAnalyzer",
    "VariantCaller", 
    "Visualizer",
    "FileHandler",
    "QualityControl"
]

if GUI_AVAILABLE:
    __all__.append("GenomicsToolkitGUI")