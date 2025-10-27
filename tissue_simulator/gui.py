"""
GUI module for interactive tissue simulation.
"""

import sys
import json
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QSlider, QPushButton, QTextEdit, QGroupBox, QRadioButton,
    QMessageBox, QFileDialog, QProgressBar, QTabWidget, QTableWidget,
    QTableWidgetItem, QCheckBox
)
from PyQt5.QtCore import Qt, QTimer, QThread, pyqtSignal
from PyQt5.QtGui import QFont
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

from .tissue import TissueSection
from .packing import SpherePacker


class PackingThread(QThread):
    """Thread for running cell packing without blocking GUI."""
    progress = pyqtSignal(int, int)  # cells_placed, total_attempts
    finished = pyqtSignal(list)  # List of cells
    
    def __init__(self, packer, max_attempts):
        super().__init__()
        self.packer = packer
        self.max_attempts = max_attempts
        
    def run(self):
        cells = self.packer.pack_with_progress(
            max_attempts=self.max_attempts,
            callback=lambda placed, attempts: self.progress.emit(placed, attempts)
        )
        self.finished.emit(cells)


class TissueViewer3D(FigureCanvas):
    """3D visualization canvas for tissue sections."""
    
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = self.fig.add_subplot(111, projection='3d')
        super().__init__(self.fig)
        self.setParent(parent)
        
        self.tissue = None
        self.elevation = 20
        self.azimuth = 45
        
    def plot_tissue(self, tissue: TissueSection, 
                   show_boundary: bool = True):
        """Plot tissue section in 3D."""
        self.tissue = tissue
        self.ax.clear()
        
        if not tissue.cells:
            self.ax.text2D(0.5, 0.5, 'No cells generated yet', 
                          transform=self.ax.transAxes,
                          ha='center', va='center', fontsize=14)
            self.draw()
            return
        
        # Color map for cell types
        cell_types = list(set(c.cell_type for c in tissue.cells))
        colors = matplotlib.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        # Plot cells (reduced resolution for performance)
        for cell in tissue.cells:
            u = np.linspace(0, 2 * np.pi, 15)
            v = np.linspace(0, np.pi, 15)
            x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
            y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
            z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
            
            color = color_map[cell.cell_type]
            alpha = 0.3 if cell.is_boundary else 0.6
            
            # Convert color to proper format for plot_surface
            self.ax.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)),
                               alpha=alpha, linewidth=0, antialiased=True, shade=False)
        
        # Draw boundary box
        if show_boundary:
            self._draw_boundary_box(tissue)
        
        # Set labels and limits
        self.ax.set_xlabel('Width (μm)')
        self.ax.set_ylabel('Height (μm)')
        self.ax.set_zlabel('Thickness (μm)')
        self.ax.set_xlim(0, tissue.width)
        self.ax.set_ylim(0, tissue.height)
        self.ax.set_zlim(0, tissue.thickness)
        
        # Set viewing angle
        self.ax.view_init(elev=self.elevation, azim=self.azimuth)
        
        # Add legend
        legend_elements = [
            matplotlib.lines.Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor=color_map[ct], 
                                   markersize=10, label=ct)
            for ct in cell_types
        ]
        self.ax.legend(handles=legend_elements, loc='upper right')
        
        self.ax.set_title(f'3D Tissue Section: {len(tissue.cells)} cells')
        self.fig.tight_layout()
        self.draw()
    
    def _draw_boundary_box(self, tissue):
        """Draw the tissue boundary box."""
        edges = [
            [[0, tissue.width], [0, 0], [0, 0]],
            [[0, tissue.width], [tissue.height, tissue.height], [0, 0]],
            [[0, tissue.width], [0, 0], [tissue.thickness, tissue.thickness]],
            [[0, tissue.width], [tissue.height, tissue.height], [tissue.thickness, tissue.thickness]],
            [[0, 0], [0, tissue.height], [0, 0]],
            [[tissue.width, tissue.width], [0, tissue.height], [0, 0]],
            [[0, 0], [0, tissue.height], [tissue.thickness, tissue.thickness]],
            [[tissue.width, tissue.width], [0, tissue.height], [tissue.thickness, tissue.thickness]],
            [[0, 0], [0, 0], [0, tissue.thickness]],
            [[tissue.width, tissue.width], [0, 0], [0, tissue.thickness]],
            [[0, 0], [tissue.height, tissue.height], [0, tissue.thickness]],
            [[tissue.width, tissue.width], [tissue.height, tissue.height], [0, tissue.thickness]],
        ]
        
        for edge in edges:
            self.ax.plot3D(*edge, 'k-', linewidth=1, alpha=0.3)
    
    def update_view_angle(self, elevation, azimuth):
        """Update 3D viewing angle."""
        self.elevation = elevation
        self.azimuth = azimuth
        if self.tissue:
            self.plot_tissue(self.tissue)


class TissueSimulatorGUI(QMainWindow):
    """Main GUI window for tissue simulator."""
    
    def __init__(self):
        super().__init__()
        self.tissue = None
        self.packing_thread = None
        self.init_ui()
        
    def init_ui(self):
        """Initialize user interface."""
        self.setWindowTitle('3D Tissue Simulator')
        self.setGeometry(100, 100, 1400, 900)
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Left panel for controls
        control_panel = self.create_control_panel()
        main_layout.addWidget(control_panel, stretch=1)
        
        # Right panel for visualization
        viz_panel = self.create_visualization_panel()
        main_layout.addWidget(viz_panel, stretch=2)
        
    def create_control_panel(self):
        """Create left control panel."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        # Tissue dimensions group
        dims_group = QGroupBox("Tissue Dimensions (μm)")
        dims_layout = QVBoxLayout()
        
        self.height_slider = self.create_slider_with_label(
            "Height:", 100, 1000, 500, dims_layout
        )
        self.width_slider = self.create_slider_with_label(
            "Width:", 100, 1000, 500, dims_layout
        )
        self.thickness_slider = self.create_slider_with_label(
            "Thickness:", 20, 200, 100, dims_layout
        )
        
        dims_group.setLayout(dims_layout)
        layout.addWidget(dims_group)
        
        # Cell radii configuration group
        radii_group = QGroupBox("Cell Radii Configuration")
        radii_layout = QVBoxLayout()
        
        # Radio buttons for input mode
        self.simple_mode_radio = QRadioButton("Simple Range")
        self.json_mode_radio = QRadioButton("JSON Cell Types")
        self.simple_mode_radio.setChecked(True)
        self.simple_mode_radio.toggled.connect(self.toggle_radii_mode)
        
        radii_layout.addWidget(self.simple_mode_radio)
        radii_layout.addWidget(self.json_mode_radio)
        
        # Simple range sliders
        self.simple_widget = QWidget()
        simple_layout = QVBoxLayout()
        self.min_radius_slider = self.create_slider_with_label(
            "Min Radius:", 1, 20, 5, simple_layout
        )
        self.max_radius_slider = self.create_slider_with_label(
            "Max Radius:", 5, 30, 10, simple_layout
        )
        self.simple_widget.setLayout(simple_layout)
        radii_layout.addWidget(self.simple_widget)
        
        # JSON text input
        self.json_widget = QWidget()
        json_layout = QVBoxLayout()
        json_layout.addWidget(QLabel("JSON Format:"))
        self.json_text = QTextEdit()
        self.json_text.setPlaceholderText(
            '{\n  "epithelial": [5, 10],\n  "stromal": [8, 15],\n  "immune": [3, 6]\n}'
        )
        self.json_text.setMaximumHeight(150)
        json_layout.addWidget(self.json_text)
        self.json_widget.setLayout(json_layout)
        self.json_widget.setVisible(False)
        radii_layout.addWidget(self.json_widget)
        
        radii_group.setLayout(radii_layout)
        layout.addWidget(radii_group)
        
        # Packing parameters group
        packing_group = QGroupBox("Packing Parameters")
        packing_layout = QVBoxLayout()
        
        self.max_attempts_slider = self.create_slider_with_label(
            "Max Attempts:", 100, 5000, 1000, packing_layout
        )
        self.min_spacing_slider = self.create_slider_with_label(
            "Min Spacing:", 0, 50, 5, packing_layout, scale=10
        )
        
        self.boundary_cells_checkbox = QCheckBox("Allow Boundary Cells")
        self.boundary_cells_checkbox.setChecked(True)
        packing_layout.addWidget(self.boundary_cells_checkbox)
        
        packing_group.setLayout(packing_layout)
        layout.addWidget(packing_group)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        self.progress_label = QLabel("")
        layout.addWidget(self.progress_label)
        
        # Action buttons
        button_layout = QHBoxLayout()
        
        self.generate_btn = QPushButton("Generate Tissue")
        self.generate_btn.clicked.connect(self.generate_tissue)
        button_layout.addWidget(self.generate_btn)
        
        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self.clear_tissue)
        button_layout.addWidget(self.clear_btn)
        
        layout.addLayout(button_layout)
        
        # Export button
        self.export_btn = QPushButton("Export to CSV")
        self.export_btn.clicked.connect(self.export_tissue)
        self.export_btn.setEnabled(False)
        layout.addWidget(self.export_btn)
        
        layout.addStretch()
        
        return panel
    
    def create_visualization_panel(self):
        """Create right visualization panel."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        # Tab widget for different views
        tabs = QTabWidget()
        
        # 3D Viewer tab
        viewer_tab = QWidget()
        viewer_layout = QVBoxLayout(viewer_tab)
        
        # 3D canvas
        self.viewer_3d = TissueViewer3D(viewer_tab, width=8, height=6)
        viewer_layout.addWidget(self.viewer_3d)
        
        # View controls
        view_controls = QWidget()
        view_layout = QHBoxLayout(view_controls)
        
        self.elevation_slider = self.create_slider_with_label(
            "Elevation:", -90, 90, 20, view_layout, horizontal=True
        )
        self.azimuth_slider = self.create_slider_with_label(
            "Azimuth:", 0, 360, 45, view_layout, horizontal=True
        )
        
        self.elevation_slider.valueChanged.connect(self.update_view)
        self.azimuth_slider.valueChanged.connect(self.update_view)
        
        viewer_layout.addWidget(view_controls)
        
        tabs.addTab(viewer_tab, "3D Viewer")
        
        # Statistics tab
        stats_tab = QWidget()
        stats_layout = QVBoxLayout(stats_tab)
        
        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setFont(QFont("Courier", 10))
        stats_layout.addWidget(self.stats_text)
        
        tabs.addTab(stats_tab, "Statistics")
        
        layout.addWidget(tabs)
        
        return panel
    
    def create_slider_with_label(self, label_text, min_val, max_val, 
                                 default_val, layout, scale=1, horizontal=False):
        """Create a slider with label and value display."""
        if horizontal:
            container = QWidget()
            container_layout = QVBoxLayout(container)
            container_layout.setContentsMargins(0, 0, 0, 0)
        else:
            container = None
            container_layout = layout
        
        label = QLabel(f"{label_text} {default_val / scale:.1f}")
        container_layout.addWidget(label)
        
        slider = QSlider(Qt.Horizontal)
        slider.setMinimum(int(min_val * scale))
        slider.setMaximum(int(max_val * scale))
        slider.setValue(int(default_val * scale))
        slider.setTickPosition(QSlider.TicksBelow)
        slider.setTickInterval(int((max_val - min_val) * scale / 10))
        
        def update_label(value):
            label.setText(f"{label_text} {value / scale:.1f}")
        
        slider.valueChanged.connect(update_label)
        container_layout.addWidget(slider)
        
        if horizontal:
            layout.addWidget(container)
        
        return slider
    
    def toggle_radii_mode(self):
        """Toggle between simple and JSON radii input modes."""
        simple_mode = self.simple_mode_radio.isChecked()
        self.simple_widget.setVisible(simple_mode)
        self.json_widget.setVisible(not simple_mode)
    
    def get_cell_radii_config(self):
        """Get cell radii configuration from UI."""
        if self.simple_mode_radio.isChecked():
            min_r = self.min_radius_slider.value()
            max_r = self.max_radius_slider.value()
            
            if min_r >= max_r:
                QMessageBox.warning(
                    self, "Invalid Range",
                    "Minimum radius must be less than maximum radius."
                )
                return None
            
            return {"default": (min_r, max_r)}
        else:
            try:
                json_str = self.json_text.toPlainText()
                radii_dict = json.loads(json_str)
                
                # Validate format
                for cell_type, radii in radii_dict.items():
                    if not isinstance(radii, (list, tuple)) or len(radii) != 2:
                        raise ValueError(
                            f"Invalid format for '{cell_type}': must be [min, max]"
                        )
                    if radii[0] >= radii[1]:
                        raise ValueError(
                            f"Invalid range for '{cell_type}': min must be < max"
                        )
                
                # Convert lists to tuples
                return {k: tuple(v) for k, v in radii_dict.items()}
                
            except json.JSONDecodeError as e:
                QMessageBox.warning(
                    self, "JSON Error",
                    f"Invalid JSON format:\n{str(e)}"
                )
                return None
            except ValueError as e:
                QMessageBox.warning(
                    self, "Configuration Error",
                    str(e)
                )
                return None
    
    def generate_tissue(self):
        """Generate tissue with current parameters."""
        # Get parameters
        height = self.height_slider.value()
        width = self.width_slider.value()
        thickness = self.thickness_slider.value()
        
        cell_radii = self.get_cell_radii_config()
        if cell_radii is None:
            return
        
        max_attempts = self.max_attempts_slider.value()
        min_spacing = self.min_spacing_slider.value() / 10.0
        allow_boundary = self.boundary_cells_checkbox.isChecked()
        
        # Create tissue
        self.tissue = TissueSection(
            height=height,
            width=width,
            thickness=thickness,
            cell_radii=cell_radii
        )
        
        # Create packer
        packer = SpherePacker(
            bounds=self.tissue.get_bounds(),
            cell_radii_config=cell_radii,
            min_spacing=min_spacing,
            allow_boundary_cells=allow_boundary
        )
        
        # Disable controls during generation
        self.generate_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.progress_label.setText("Generating cells...")
        
        # Start packing thread
        self.packing_thread = PackingThread(packer, max_attempts)
        self.packing_thread.progress.connect(self.update_progress)
        self.packing_thread.finished.connect(self.packing_finished)
        self.packing_thread.start()
    
    def update_progress(self, cells_placed, total_attempts):
        """Update progress bar during cell generation."""
        progress = min(100, int(cells_placed / 10))  # Rough estimate
        self.progress_bar.setValue(progress)
        self.progress_label.setText(
            f"Cells placed: {cells_placed} (Attempts: {total_attempts})"
        )
    
    def packing_finished(self, cells):
        """Handle completion of cell packing."""
        self.tissue.cells = cells
        
        # Update UI
        self.progress_bar.setVisible(False)
        self.progress_label.setText(f"Generation complete: {len(cells)} cells")
        self.generate_btn.setEnabled(True)
        self.export_btn.setEnabled(True)
        
        # Update visualization
        self.viewer_3d.plot_tissue(self.tissue)
        
        # Update statistics
        self.update_statistics()
    
    def update_view(self):
        """Update 3D view angle."""
        elevation = self.elevation_slider.value()
        azimuth = self.azimuth_slider.value()
        self.viewer_3d.update_view_angle(elevation, azimuth)
    
    def update_statistics(self):
        """Update statistics display."""
        if self.tissue is None or not self.tissue.cells:
            self.stats_text.setPlainText("No tissue data available.")
            return
        
        stats = self.tissue.get_cell_statistics()
        
        text = "=== Tissue Statistics ===\n\n"
        text += f"Tissue Dimensions: {self.tissue.width} x {self.tissue.height} x {self.tissue.thickness} μm\n"
        text += f"Tissue Volume: {self.tissue.width * self.tissue.height * self.tissue.thickness:.0f} μm³\n\n"
        
        text += f"Total Cells: {stats['total_cells']}\n"
        text += f"Interior Cells: {stats['interior_cells']}\n"
        text += f"Boundary Cells: {stats['boundary_cells']}\n"
        text += f"Packing Fraction: {stats['packing_fraction']:.3f}\n\n"
        
        text += "=== Cell Types ===\n\n"
        for cell_type, count in stats['cell_types'].items():
            avg_radius = stats['avg_radii'][cell_type]
            text += f"{cell_type}:\n"
            text += f"  Count: {count}\n"
            text += f"  Avg Radius: {avg_radius:.2f} μm\n\n"
        
        self.stats_text.setPlainText(text)
    
    def clear_tissue(self):
        """Clear current tissue."""
        if self.tissue:
            self.tissue.clear_cells()
        
        self.tissue = None
        self.viewer_3d.ax.clear()
        self.viewer_3d.draw()
        self.stats_text.clear()
        self.progress_label.setText("")
        self.export_btn.setEnabled(False)
    
    def export_tissue(self):
        """Export tissue data to CSV."""
        if self.tissue is None or not self.tissue.cells:
            QMessageBox.warning(
                self, "No Data",
                "No tissue data to export."
            )
            return
        
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Tissue Data", "", "CSV Files (*.csv)"
        )
        
        if filename:
            try:
                self.tissue.export_to_csv(filename)
                QMessageBox.information(
                    self, "Export Successful",
                    f"Tissue data exported to:\n{filename}"
                )
            except Exception as e:
                QMessageBox.critical(
                    self, "Export Failed",
                    f"Failed to export data:\n{str(e)}"
                )


def main():
    """Launch the GUI application."""
    app = QApplication(sys.argv)
    app.setStyle('Fusion')  # Modern look
    
    window = TissueSimulatorGUI()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
