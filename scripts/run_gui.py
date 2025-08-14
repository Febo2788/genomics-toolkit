#!/usr/bin/env python3
"""
Launch script for GenomicsToolkit GUI
"""

import sys
import os
from pathlib import Path

# Add the parent directory to the path so we can import genomics_toolkit
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from genomics_toolkit.gui import main
    
    if __name__ == "__main__":
        main()
        
except ImportError as e:
    print(f"Error importing GenomicsToolkit GUI: {e}")
    print("Please install the required dependencies:")
    print("pip install PyQt6 pyqtgraph matplotlib")
    sys.exit(1)
except Exception as e:
    print(f"Error launching GUI: {e}")
    sys.exit(1)