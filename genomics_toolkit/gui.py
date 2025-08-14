"""
PyQt6 GUI Interface for GenomicsToolkit

Modern graphical interface for sequence analysis, variant calling, and visualization.
"""

import sys
import os
import json
from pathlib import Path
from typing import Optional, List, Dict, Any
import traceback

import numpy as np
import pandas as pd

# Set matplotlib backend before importing pyplot
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# Try to import Qt5 backend, fallback to basic canvas
try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
except ImportError:
    # Fallback: create a simple widget that can display matplotlib figures
    from PyQt6.QtWidgets import QLabel
    from PyQt6.QtGui import QPixmap
    from PyQt6.QtCore import QByteArray
    import io
    
    class FigureCanvas(QLabel):
        def __init__(self, figure):
            super().__init__()
            self.figure = figure
            self.draw()
            
        def draw(self):
            buf = io.BytesIO()
            self.figure.savefig(buf, format='png', dpi=100, bbox_inches='tight')
            buf.seek(0)
            pixmap = QPixmap()
            pixmap.loadFromData(buf.getvalue())
            self.setPixmap(pixmap)
            buf.close()

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QGridLayout, QTabWidget, QLabel, QPushButton, QTextEdit, QLineEdit,
    QFileDialog, QProgressBar, QTableWidget, QTableWidgetItem, QGroupBox,
    QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox, QMessageBox, QSplitter,
    QScrollArea, QFrame, QSizePolicy, QToolBar, QStatusBar, QMenuBar
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QSize
from PyQt6.QtGui import QPixmap, QFont, QIcon, QAction, QPalette, QColor

from .sequence_analysis import SequenceAnalyzer
from .variant_calling import VariantCaller
from .visualization import Visualizer
from .utils import FileHandler, QualityControl


class AnalysisWorker(QThread):
    """Background thread for running analysis tasks."""
    
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    
    def __init__(self, analysis_type: str, file_path: str, parameters: dict):
        super().__init__()
        self.analysis_type = analysis_type
        self.file_path = file_path
        self.parameters = parameters
        
    def run(self):
        """Run the analysis in background thread."""
        try:
            analyzer = SequenceAnalyzer()
            caller = VariantCaller()
            file_handler = FileHandler()
            
            if self.analysis_type == "sequence_analysis":
                self.status.emit("Loading sequences...")
                self.progress.emit(10)
                
                sequences = list(analyzer.parse_fasta(self.file_path))
                self.progress.emit(30)
                
                results = {"sequences": [], "stats_summary": {}}
                
                for i, (header, seq) in enumerate(sequences):
                    self.status.emit(f"Analyzing sequence {i+1}/{len(sequences)}")
                    
                    # Basic stats
                    stats = analyzer.basic_stats(seq)
                    
                    # ORF detection
                    orfs = analyzer.find_orfs(seq, self.parameters.get("min_orf_length", 100))
                    
                    # Composition analysis
                    composition = analyzer.composition_analysis(seq, self.parameters.get("window_size", 100))
                    
                    seq_result = {
                        "header": header,
                        "length": len(seq),
                        "stats": stats,
                        "orfs": orfs,
                        "composition": composition,
                        "sequence": seq[:1000] if len(seq) > 1000 else seq  # Truncate for display
                    }
                    results["sequences"].append(seq_result)
                    
                    progress = 30 + (i + 1) / len(sequences) * 60
                    self.progress.emit(int(progress))
                
                # Summary statistics
                all_lengths = [len(seq) for _, seq in sequences]
                all_gc = [analyzer.basic_stats(seq)["gc_content"] for _, seq in sequences]
                
                results["stats_summary"] = {
                    "total_sequences": len(sequences),
                    "total_length": sum(all_lengths),
                    "mean_length": np.mean(all_lengths),
                    "n50": analyzer.calculate_n50([seq for _, seq in sequences]),
                    "mean_gc": np.mean(all_gc)
                }
                
                self.progress.emit(100)
                self.status.emit("Analysis complete!")
                self.finished.emit(results)
                
            elif self.analysis_type == "variant_calling":
                self.status.emit("Loading reference sequence...")
                self.progress.emit(10)
                
                sequences = list(analyzer.parse_fasta(self.file_path))
                if not sequences:
                    raise ValueError("No sequences found in file")
                
                ref_seq = sequences[0][1]
                self.progress.emit(30)
                
                self.status.emit("Simulating variants...")
                variants = caller.simulate_variants(
                    ref_seq, 
                    self.parameters.get("num_snps", 50),
                    self.parameters.get("num_indels", 10)
                )
                self.progress.emit(60)
                
                self.status.emit("Filtering variants...")
                filtered = caller.filter_variants(
                    min_quality=self.parameters.get("min_quality", 20),
                    min_allele_freq=self.parameters.get("min_allele_freq", 0.05)
                )
                self.progress.emit(80)
                
                self.status.emit("Annotating variants...")
                annotated = caller.annotate_variants(filtered)
                
                # Calculate statistics
                pop_stats = caller.calculate_population_stats(annotated)
                
                results = {
                    "variants": annotated,
                    "statistics": pop_stats,
                    "reference_length": len(ref_seq)
                }
                
                self.progress.emit(100)
                self.status.emit("Variant calling complete!")
                self.finished.emit(results)
                
        except Exception as e:
            self.error.emit(f"Analysis failed: {str(e)}")


class PlotWidget(QWidget):
    """Custom widget for displaying matplotlib plots."""
    
    def __init__(self):
        super().__init__()
        self.figure = Figure(figsize=(10, 6))
        self.canvas = FigureCanvas(self.figure)
        
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        
    def clear(self):
        """Clear the plot."""
        self.figure.clear()
        self.canvas.draw()
        
    def plot_composition(self, composition_data):
        """Plot sequence composition."""
        self.figure.clear()
        ax1 = self.figure.add_subplot(211)
        ax2 = self.figure.add_subplot(212)
        
        positions = composition_data["positions"]
        
        # GC content
        ax1.plot(positions, composition_data["gc_content"], 'b-', linewidth=2)
        ax1.fill_between(positions, composition_data["gc_content"], alpha=0.3)
        ax1.set_ylabel('GC Content (%)')
        ax1.set_title('Sequence Composition')
        ax1.grid(True, alpha=0.3)
        
        # AT content
        ax2.plot(positions, composition_data["at_content"], 'r-', linewidth=2)
        ax2.fill_between(positions, composition_data["at_content"], alpha=0.3)
        ax2.set_ylabel('AT Content (%)')
        ax2.set_xlabel('Position (bp)')
        ax2.grid(True, alpha=0.3)
        
        self.figure.tight_layout()
        self.canvas.draw()
        
    def plot_variants(self, variants):
        """Plot variant statistics."""
        if not variants:
            return
            
        self.figure.clear()
        
        # Create 2x2 subplot
        ax1 = self.figure.add_subplot(221)
        ax2 = self.figure.add_subplot(222)
        ax3 = self.figure.add_subplot(223)
        ax4 = self.figure.add_subplot(224)
        
        # Variant types
        types = []
        for var in variants:
            if var.is_snp():
                types.append('SNP')
            elif var.is_insertion():
                types.append('Insertion')
            elif var.is_deletion():
                types.append('Deletion')
        
        type_counts = pd.Series(types).value_counts()
        ax1.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%')
        ax1.set_title('Variant Types')
        
        # Quality distribution
        qualities = [var.quality for var in variants]
        ax2.hist(qualities, bins=20, alpha=0.7, color='skyblue')
        ax2.set_xlabel('Quality Score')
        ax2.set_ylabel('Count')
        ax2.set_title('Quality Distribution')
        
        # Allele frequencies
        allele_freqs = [var.info.get('AF', 0) for var in variants]
        ax3.hist(allele_freqs, bins=20, alpha=0.7, color='lightgreen')
        ax3.set_xlabel('Allele Frequency')
        ax3.set_ylabel('Count')
        ax3.set_title('Allele Frequency Distribution')
        
        # Position vs Quality
        positions = [var.pos for var in variants]
        ax4.scatter(positions, qualities, alpha=0.6)
        ax4.set_xlabel('Position')
        ax4.set_ylabel('Quality')
        ax4.set_title('Quality vs Position')
        
        self.figure.tight_layout()
        self.canvas.draw()


class GenomicsToolkitGUI(QMainWindow):
    """Main GUI window for GenomicsToolkit."""
    
    def __init__(self):
        super().__init__()
        self.current_results = None
        self.init_ui()
        
    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle("GenomicsToolkit - Bioinformatics Analysis Suite")
        self.setGeometry(100, 100, 1400, 900)
        
        # Set application style
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f5f5f5;
            }
            QTabWidget::pane {
                border: 1px solid #c0c0c0;
                background-color: white;
            }
            QTabBar::tab {
                background-color: #e0e0e0;
                padding: 8px 16px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: white;
                border-bottom: 2px solid #0078d4;
            }
            QPushButton {
                background-color: #0078d4;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #106ebe;
            }
            QPushButton:pressed {
                background-color: #005a9e;
            }
            QGroupBox {
                font-weight: bold;
                border: 2px solid #cccccc;
                border-radius: 5px;
                margin-top: 1ex;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        
        # Create menu bar
        self.create_menu_bar()
        
        # Create status bar
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready")
        
        # Create main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        
        layout = QVBoxLayout(main_widget)
        
        # Create tab widget
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # Add tabs
        self.create_input_tab()
        self.create_results_tab()
        self.create_visualization_tab()
        self.create_settings_tab()
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
    def create_menu_bar(self):
        """Create the menu bar."""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('File')
        
        open_action = QAction('Open FASTA/FASTQ', self)
        open_action.triggered.connect(self.select_file)
        file_menu.addAction(open_action)
        
        save_action = QAction('Save Results', self)
        save_action.triggered.connect(self.save_results)
        file_menu.addAction(save_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Help menu
        help_menu = menubar.addMenu('Help')
        
        about_action = QAction('About', self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
    def create_input_tab(self):
        """Create the input/analysis tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # File selection group
        file_group = QGroupBox("Input File")
        file_layout = QHBoxLayout(file_group)
        
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select a FASTA or FASTQ file...")
        file_layout.addWidget(self.file_path_edit)
        
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.select_file)
        file_layout.addWidget(browse_btn)
        
        layout.addWidget(file_group)
        
        # Analysis type group
        analysis_group = QGroupBox("Analysis Type")
        analysis_layout = QVBoxLayout(analysis_group)
        
        self.analysis_combo = QComboBox()
        self.analysis_combo.addItems(["Sequence Analysis", "Variant Calling"])
        self.analysis_combo.currentTextChanged.connect(self.on_analysis_type_changed)
        analysis_layout.addWidget(self.analysis_combo)
        
        layout.addWidget(analysis_group)
        
        # Parameters group
        self.params_group = QGroupBox("Parameters")
        self.params_layout = QGridLayout(self.params_group)
        self.create_sequence_params()
        layout.addWidget(self.params_group)
        
        # Run button
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        self.run_btn.setFixedHeight(40)
        layout.addWidget(self.run_btn)
        
        layout.addStretch()
        self.tab_widget.addTab(tab, "Analysis")
        
    def create_sequence_params(self):
        """Create parameter controls for sequence analysis."""
        # Clear existing widgets
        for i in reversed(range(self.params_layout.count())):
            self.params_layout.itemAt(i).widget().setParent(None)
        
        if self.analysis_combo.currentText() == "Sequence Analysis":
            self.params_layout.addWidget(QLabel("Minimum ORF Length:"), 0, 0)
            self.min_orf_spin = QSpinBox()
            self.min_orf_spin.setRange(30, 1000)
            self.min_orf_spin.setValue(100)
            self.min_orf_spin.setSuffix(" bp")
            self.params_layout.addWidget(self.min_orf_spin, 0, 1)
            
            self.params_layout.addWidget(QLabel("Window Size:"), 1, 0)
            self.window_size_spin = QSpinBox()
            self.window_size_spin.setRange(50, 500)
            self.window_size_spin.setValue(100)
            self.window_size_spin.setSuffix(" bp")
            self.params_layout.addWidget(self.window_size_spin, 1, 1)
            
        else:  # Variant Calling
            self.params_layout.addWidget(QLabel("Number of SNPs:"), 0, 0)
            self.num_snps_spin = QSpinBox()
            self.num_snps_spin.setRange(10, 500)
            self.num_snps_spin.setValue(50)
            self.params_layout.addWidget(self.num_snps_spin, 0, 1)
            
            self.params_layout.addWidget(QLabel("Number of Indels:"), 1, 0)
            self.num_indels_spin = QSpinBox()
            self.num_indels_spin.setRange(5, 100)
            self.num_indels_spin.setValue(10)
            self.params_layout.addWidget(self.num_indels_spin, 1, 1)
            
            self.params_layout.addWidget(QLabel("Min Quality:"), 2, 0)
            self.min_quality_spin = QDoubleSpinBox()
            self.min_quality_spin.setRange(10.0, 60.0)
            self.min_quality_spin.setValue(20.0)
            self.params_layout.addWidget(self.min_quality_spin, 2, 1)
            
            self.params_layout.addWidget(QLabel("Min Allele Freq:"), 3, 0)
            self.min_af_spin = QDoubleSpinBox()
            self.min_af_spin.setRange(0.01, 0.5)
            self.min_af_spin.setValue(0.05)
            self.min_af_spin.setDecimals(3)
            self.params_layout.addWidget(self.min_af_spin, 3, 1)
    
    def create_results_tab(self):
        """Create the results display tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Results summary
        self.summary_text = QTextEdit()
        self.summary_text.setMaximumHeight(200)
        self.summary_text.setFont(QFont("Consolas", 9))
        layout.addWidget(QLabel("Analysis Summary:"))
        layout.addWidget(self.summary_text)
        
        # Results table
        self.results_table = QTableWidget()
        layout.addWidget(QLabel("Detailed Results:"))
        layout.addWidget(self.results_table)
        
        # Export button
        export_btn = QPushButton("Export Results")
        export_btn.clicked.connect(self.export_results)
        layout.addWidget(export_btn)
        
        self.tab_widget.addTab(tab, "Results")
        
    def create_visualization_tab(self):
        """Create the visualization tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Plot controls
        controls_layout = QHBoxLayout()
        
        plot_type_combo = QComboBox()
        plot_type_combo.addItems(["Sequence Composition", "Variant Statistics", "ORF Distribution"])
        plot_type_combo.currentTextChanged.connect(self.update_plot)
        controls_layout.addWidget(QLabel("Plot Type:"))
        controls_layout.addWidget(plot_type_combo)
        
        save_plot_btn = QPushButton("Save Plot")
        save_plot_btn.clicked.connect(self.save_plot)
        controls_layout.addWidget(save_plot_btn)
        
        controls_layout.addStretch()
        layout.addLayout(controls_layout)
        
        # Plot widget
        self.plot_widget = PlotWidget()
        layout.addWidget(self.plot_widget)
        
        self.tab_widget.addTab(tab, "Visualization")
        
    def create_settings_tab(self):
        """Create the settings tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        settings_group = QGroupBox("Application Settings")
        settings_layout = QGridLayout(settings_group)
        
        settings_layout.addWidget(QLabel("Output Directory:"), 0, 0)
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setText(str(Path.cwd() / "output"))
        settings_layout.addWidget(self.output_dir_edit, 0, 1)
        
        browse_output_btn = QPushButton("Browse")
        browse_output_btn.clicked.connect(self.select_output_dir)
        settings_layout.addWidget(browse_output_btn, 0, 2)
        
        settings_layout.addWidget(QLabel("Max Workers:"), 1, 0)
        self.max_workers_spin = QSpinBox()
        self.max_workers_spin.setRange(1, 16)
        self.max_workers_spin.setValue(4)
        settings_layout.addWidget(self.max_workers_spin, 1, 1)
        
        layout.addWidget(settings_group)
        layout.addStretch()
        
        self.tab_widget.addTab(tab, "Settings")
    
    def on_analysis_type_changed(self):
        """Handle analysis type change."""
        self.create_sequence_params()
        
    def select_file(self):
        """Open file dialog to select input file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select FASTA or FASTQ file",
            "",
            "FASTA files (*.fasta *.fa *.fas);;FASTQ files (*.fastq *.fq);;All files (*.*)"
        )
        
        if file_path:
            self.file_path_edit.setText(file_path)
            
            # Validate file
            file_handler = FileHandler()
            validation = file_handler.validate_file_integrity(file_path)
            
            if validation["valid"]:
                self.status_bar.showMessage(f"File loaded: {Path(file_path).name} ({validation['file_format']})")
            else:
                QMessageBox.warning(self, "File Validation", f"File validation failed: {validation['error']}")
    
    def select_output_dir(self):
        """Select output directory."""
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir_edit.setText(dir_path)
    
    def run_analysis(self):
        """Run the selected analysis."""
        file_path = self.file_path_edit.text()
        
        if not file_path or not Path(file_path).exists():
            QMessageBox.warning(self, "Error", "Please select a valid input file.")
            return
        
        # Get parameters
        if self.analysis_combo.currentText() == "Sequence Analysis":
            parameters = {
                "min_orf_length": self.min_orf_spin.value(),
                "window_size": self.window_size_spin.value()
            }
            analysis_type = "sequence_analysis"
        else:
            parameters = {
                "num_snps": self.num_snps_spin.value(),
                "num_indels": self.num_indels_spin.value(),
                "min_quality": self.min_quality_spin.value(),
                "min_allele_freq": self.min_af_spin.value()
            }
            analysis_type = "variant_calling"
        
        # Disable UI during analysis
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        
        # Start analysis in background thread
        self.worker = AnalysisWorker(analysis_type, file_path, parameters)
        self.worker.progress.connect(self.progress_bar.setValue)
        self.worker.status.connect(self.status_bar.showMessage)
        self.worker.finished.connect(self.on_analysis_finished)
        self.worker.error.connect(self.on_analysis_error)
        self.worker.start()
    
    def on_analysis_finished(self, results):
        """Handle analysis completion."""
        self.current_results = results
        self.run_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        # Update results display
        self.update_results_display()
        self.update_plot()
        
        # Switch to results tab
        self.tab_widget.setCurrentIndex(1)
        
        QMessageBox.information(self, "Analysis Complete", "Analysis completed successfully!")
    
    def on_analysis_error(self, error_msg):
        """Handle analysis error."""
        self.run_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.status_bar.showMessage("Analysis failed")
        
        QMessageBox.critical(self, "Analysis Error", error_msg)
    
    def update_results_display(self):
        """Update the results display."""
        if not self.current_results:
            return
        
        # Update summary
        if "stats_summary" in self.current_results:
            summary = self.current_results["stats_summary"]
            summary_text = f"""
SEQUENCE ANALYSIS SUMMARY
========================
Total sequences: {summary['total_sequences']}
Total length: {summary['total_length']:,} bp
Mean length: {summary['mean_length']:.0f} bp
N50: {summary['n50']:,} bp
Mean GC content: {summary['mean_gc']:.1f}%
"""
            self.summary_text.setText(summary_text)
            
            # Update table with sequence details
            sequences = self.current_results.get("sequences", [])
            self.results_table.setRowCount(len(sequences))
            self.results_table.setColumnCount(5)
            self.results_table.setHorizontalHeaderLabels([
                "Sequence", "Length", "GC%", "ORFs", "AT%"
            ])
            
            for i, seq_data in enumerate(sequences):
                self.results_table.setItem(i, 0, QTableWidgetItem(seq_data["header"][:50]))
                self.results_table.setItem(i, 1, QTableWidgetItem(f"{seq_data['length']:,}"))
                self.results_table.setItem(i, 2, QTableWidgetItem(f"{seq_data['stats']['gc_content']:.1f}"))
                self.results_table.setItem(i, 3, QTableWidgetItem(str(len(seq_data['orfs']))))
                self.results_table.setItem(i, 4, QTableWidgetItem(f"{seq_data['stats']['at_content']:.1f}"))
                
        elif "statistics" in self.current_results:
            stats = self.current_results["statistics"]
            summary_text = f"""
VARIANT CALLING SUMMARY
======================
Total variants: {stats['total_variants']}
SNPs: {stats['snp_count']}
Indels: {stats['indel_count']}
SNP ratio: {stats['snp_ratio']:.2f}
Ti/Tv ratio: {stats['ti_tv_ratio']:.2f}
Mean allele frequency: {stats['mean_allele_freq']:.3f}
"""
            self.summary_text.setText(summary_text)
            
            # Update table with variant details
            variants = self.current_results.get("variants", [])
            self.results_table.setRowCount(min(len(variants), 100))  # Limit to 100 for display
            self.results_table.setColumnCount(6)
            self.results_table.setHorizontalHeaderLabels([
                "Position", "Type", "Ref", "Alt", "Quality", "Allele Freq"
            ])
            
            for i, variant in enumerate(variants[:100]):
                var_type = "SNP" if variant.is_snp() else "INDEL"
                self.results_table.setItem(i, 0, QTableWidgetItem(str(variant.pos)))
                self.results_table.setItem(i, 1, QTableWidgetItem(var_type))
                self.results_table.setItem(i, 2, QTableWidgetItem(variant.ref))
                self.results_table.setItem(i, 3, QTableWidgetItem(variant.alt))
                self.results_table.setItem(i, 4, QTableWidgetItem(f"{variant.quality:.1f}"))
                self.results_table.setItem(i, 5, QTableWidgetItem(f"{variant.info.get('AF', 0):.3f}"))
        
        self.results_table.resizeColumnsToContents()
    
    def update_plot(self):
        """Update the visualization plot."""
        if not self.current_results:
            return
        
        plot_type = self.sender().currentText() if hasattr(self.sender(), 'currentText') else "Sequence Composition"
        
        if plot_type == "Sequence Composition" and "sequences" in self.current_results:
            sequences = self.current_results["sequences"]
            if sequences:
                composition = sequences[0]["composition"]  # Plot first sequence
                self.plot_widget.plot_composition(composition)
                
        elif plot_type == "Variant Statistics" and "variants" in self.current_results:
            variants = self.current_results["variants"]
            self.plot_widget.plot_variants(variants)
    
    def save_plot(self):
        """Save the current plot."""
        if not self.current_results:
            QMessageBox.warning(self, "No Data", "No plot data available to save.")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Save Plot",
            "plot.png",
            "PNG files (*.png);;PDF files (*.pdf);;All files (*.*)"
        )
        
        if file_path:
            self.plot_widget.figure.savefig(file_path, dpi=300, bbox_inches='tight')
            QMessageBox.information(self, "Success", f"Plot saved to {file_path}")
    
    def export_results(self):
        """Export analysis results."""
        if not self.current_results:
            QMessageBox.warning(self, "No Data", "No results available to export.")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Results",
            "results.json",
            "JSON files (*.json);;CSV files (*.csv);;All files (*.*)"
        )
        
        if file_path:
            try:
                if file_path.endswith('.json'):
                    # Convert results to JSON-serializable format
                    export_data = {}
                    for key, value in self.current_results.items():
                        if key == "variants":
                            export_data[key] = [
                                {
                                    "chrom": v.chrom,
                                    "pos": v.pos,
                                    "ref": v.ref,
                                    "alt": v.alt,
                                    "quality": v.quality,
                                    "info": v.info
                                } for v in value
                            ]
                        else:
                            export_data[key] = value
                    
                    with open(file_path, 'w') as f:
                        json.dump(export_data, f, indent=2, default=str)
                        
                elif file_path.endswith('.csv'):
                    # Export table data as CSV
                    table_data = []
                    for row in range(self.results_table.rowCount()):
                        row_data = []
                        for col in range(self.results_table.columnCount()):
                            item = self.results_table.item(row, col)
                            row_data.append(item.text() if item else "")
                        table_data.append(row_data)
                    
                    df = pd.DataFrame(table_data, columns=[
                        self.results_table.horizontalHeaderItem(i).text() 
                        for i in range(self.results_table.columnCount())
                    ])
                    df.to_csv(file_path, index=False)
                
                QMessageBox.information(self, "Success", f"Results exported to {file_path}")
                
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to export results: {str(e)}")
    
    def save_results(self):
        """Save results from menu."""
        self.export_results()
    
    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(self, "About GenomicsToolkit", 
                         "GenomicsToolkit v1.0.0\n\n"
                         "A comprehensive bioinformatics pipeline for sequence analysis and variant calling.\n\n"
                         "Created by Felix Borrego\n"
                         "Built with PyQt6 and Python")


def main():
    """Main function to run the GUI application."""
    app = QApplication(sys.argv)
    app.setApplicationName("GenomicsToolkit")
    app.setApplicationVersion("1.0.0")
    
    # Set application style
    app.setStyle('Fusion')
    
    # Create and show main window
    window = GenomicsToolkitGUI()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()