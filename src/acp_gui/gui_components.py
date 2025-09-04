#!/usr/bin/env python3
"""ACP GUI Components - separated from the main application logic."""

import os
from pathlib import Path
from typing import List, Optional, Dict, Any

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                              QLineEdit, QSpacerItem, QSizePolicy, QPushButton, QComboBox, QGridLayout, QDialog, 
                              QDialogButtonBox, QCheckBox)
from PySide6.QtGui import QDoubleValidator
from PySide6.QtCore import QRect, QLocale

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

# Import dependencies (these should be available in the package)
try:
    from calc.core_functionality import convolute
    from graph.plots import plot_colorbar3, add_second_axis
    from utils.helpers import random_colors
    # Note: ColorButton is defined below
except ImportError as e:
    print(f"Warning: Some dependencies not available: {e}")

HOME = Path.home()
WORKDIR = Path.cwd()




class Fragments:
    """Class to handle molecular fragments."""

    def __init__(self, indices: List[List[int]]) -> None:
        self.indices = indices
        self.colors: Optional[List[str]] = None


class FragmentDialog(QDialog):
    """Dialog for managing molecular fragments and their colors."""
    
    def __init__(self, parent: Optional[QWidget] = None, 
                 nfrag: int = 2, 
                 colors: Optional[List[str]] = None) -> None:
        super().__init__(parent)

        self.setWindowTitle("Fragment Manager")
        self.vlay = QVBoxLayout()
        grid = QGridLayout()
        self.vlay.addLayout(grid)
        
        self._nfrag = nfrag
        self._colors = colors if colors else random_colors(nfrag)
        self._cbott: List[ColorButton] = []
        self._fline: List[QLineEdit] = []

        for i in range(self._nfrag):
            message = QLabel(f"Frag.{i+1}")
            tbott = ColorButton(color=self._colors[i])
            editline = QLineEdit()
            valid = QPushButton('Set', self)
            editline.setText("prova")
            
            grid.addWidget(message, i, 0)
            grid.addWidget(tbott, i, 1)
            grid.addWidget(editline, i, 2)
            grid.addWidget(valid, i, 3)

            self._cbott.append(tbott)
            self._fline.append(editline)  

        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel 
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.accepted.connect(self._getcolors)
        self.buttonBox.rejected.connect(self.reject)

        self.vlay.addWidget(self.buttonBox)
        self.setLayout(self.vlay)

    def _getcolors(self) -> None:
        """Extract colors from color buttons."""
        self._colors = [i.color() for i in self._cbott]


class PlotCanvas(FigureCanvas):
    """Canvas for matplotlib plotting."""
    
    def __init__(self, parent: Optional[QWidget] = None, 
                 width: int = 10, height: int = 8, dpi: int = 100) -> None:
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def plot(self, data: Any, start: float, end: float, fwhm: float, spectype: str,
             colors: List[str]) -> None:
        """Plot the spectral data with given parameters."""
        try:
            convoluted = convolute(data, start, end, fwhm, spectype)
            if self.figure.axes:
                self.figure.clf()
            plot_colorbar3(self.figure, convoluted, colors=colors)
            add_second_axis(self.figure.axes[0], data[0].etenergies, data[0].obs)
            self.draw()
        except Exception as e:
            print(f"Error plotting data: {e}")


class WidgetPlot(QWidget):
    """Widget containing the plot canvas and toolbar."""
    
    def __init__(self, *args, **kwargs) -> None:
        QWidget.__init__(self, *args, **kwargs)
        self.setLayout(QVBoxLayout())
        self.canvas = PlotCanvas(self, width=10, height=8)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout().addWidget(self.toolbar)
        self.layout().addWidget(self.canvas)

    def plot(self, data: Any, start: float,
             end: float, fwhm: float,
             spectype: str, 
             colors: List[str]) -> None:
        """Plot data using the canvas."""
        self.canvas.plot(data, start,
                         end, fwhm,
                         spectype,
                         colors)


class SpectralParametersWidget(QWidget):
    """Widget for spectral parameters input."""
    
    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setup_ui()
        
    def setup_ui(self) -> None:
        """Setup the spectral parameters UI."""
        layout = QHBoxLayout()
        
        # Start frequency
        start_label = QLabel('Start:', self)
        self.start = QLineEdit(self)
        start_valid = QDoubleValidator(0.0, 1450.0, 4)
        start_valid.setLocale(QLocale('English'))
        self.start.setValidator(start_valid)
        self.start.setText('950.0')
        
        # End frequency
        end_label = QLabel('End:', self)
        self.end = QLineEdit(self)
        end_valid = QDoubleValidator(1450.0, 8000.0, 4)
        end_valid.setLocale(QLocale('English'))
        self.end.setValidator(end_valid)
        self.end.setText('1500.0')
        
        # HWHM
        hwhm_label = QLabel('HWHM:', self)
        self.hwhm = QLineEdit(self)
        self.hwhm.setValidator(QDoubleValidator(0.5, 50.0, 2))
        self.hwhm.setText('5.0')
        
        layout.addWidget(start_label)
        layout.addWidget(self.start)
        layout.addWidget(end_label)
        layout.addWidget(self.end)
        layout.addWidget(hwhm_label)
        layout.addWidget(self.hwhm)
        
        self.setLayout(layout)
        
        # Connect validators
        self.start.editingFinished.connect(self._set_min)
        self.end.editingFinished.connect(self._set_max)
    
    def _set_min(self) -> None:
        """Update minimum validator for end field."""
        minimus = float(self.start.text())
        valid = QDoubleValidator(minimus + 50., 8000., 4)
        valid.setLocale(QLocale('English'))
        self.end.setValidator(valid)

    def _set_max(self) -> None:
        """Update maximum validator for start field."""
        maxim = float(self.end.text())
        valid = QDoubleValidator(0., maxim + 50., 4)
        valid.setLocale(QLocale('English'))
        self.start.setValidator(valid)
    
    def get_parameters(self) -> Dict[str, float]:
        """Get the current spectral parameters."""
        return {
            'start': float(self.start.text()),
            'end': float(self.end.text()),
            'hwhm': float(self.hwhm.text())
        }


class PropertySelector(QComboBox):
    """ComboBox for selecting spectral properties."""
    
    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setup_properties()
    
    def setup_properties(self) -> None:
        """Setup the property options."""
        self.setGeometry(QRect(40, 40, 491, 31))
        self.addItem('')
        self.addItem('IR')
        self.addItem('VCD')
        self.addItem('RAMAN')
        self.addItem('ROA')
        
        # Initially disable all except the first
        for i in range(1, 5):
            self.model().item(i).setEnabled(False)
    
    def enable_property(self, property_name: str) -> None:
        """Enable a specific property."""
        property_map = {'IR': 1, 'VCD': 2, 'RAMAN': 3, 'ROA': 4}
        if property_name in property_map:
            self.model().item(property_map[property_name]).setEnabled(True)


class ControlButtonsWidget(QWidget):
    """Widget containing control buttons."""
    
    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setup_ui()
    
    def setup_ui(self) -> None:
        """Setup the control buttons."""
        layout = QHBoxLayout()
        
        self.groups_button = QPushButton('Load Groups', self)
        self.define_groups_button = QPushButton('Define Groups', self)
        self.interactive_fragments_button = QPushButton('Pick Fragments', self)
        self.acp_colors_cb = QCheckBox("ACP colors", self)
        self.acp_button = QPushButton('ACP Calc.', self)
        
        layout.addWidget(self.groups_button)
        layout.addWidget(self.define_groups_button)
        layout.addWidget(self.interactive_fragments_button)
        layout.addWidget(self.acp_colors_cb)
        layout.addItem(QSpacerItem(100, 10, QSizePolicy.Expanding))
        
        # Initially disable some buttons
        self.define_groups_button.setEnabled(True)
        self.interactive_fragments_button.setEnabled(True)
        self.acp_button.setEnabled(False)
        
        self.setLayout(layout)

class InteractiveFragmentDialog(QDialog):
    """Dialog for defining molecular fragments by picking atoms in Jmol."""
    
    def __init__(self, parent: Optional[QWidget] = None, 
                 jmol: Optional[Any] = None, 
                 natoms: int = 0,
                 existing_fragments: Optional[List[List[int]]] = None) -> None:
        super().__init__(parent)
        
        self.jmol = jmol
        self.natoms = natoms
        self.selected_atoms: List[int] = []  # 0-based atom indices
        self.saved_orientation: Optional[str] = None
        self.saved_selection_halos: Optional[str] = None  # Store original selectionHalos setting
        self.fragments: List[List[int]] = existing_fragments or []
        
        from utils.helpers import random_colors
        self.colors: List[str] = random_colors(max(2, len(self.fragments))) if self.fragments else random_colors(2)
        
        self.setup_ui()
        self.setup_jmol()
        
        # Initialize display
        self.update_display()
        
    def setup_ui(self) -> None:
        """Setup the dialog UI."""
        self.setWindowTitle("Define Fragments by Picking Atoms")
        self.setMinimumSize(400, 300)
        
        layout = QVBoxLayout()
        
        # Instructions
        instructions = QLabel(
            "1. Click atoms in Jmol to select them for Fragment 1\n"
            "2. Press 'Read Selection' to capture your selection\n"
            "3. Selected atoms will be highlighted with halos\n"
            "4. Fragment 2 will contain all remaining atoms"
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)
        
        # Selected atoms display
        self.selected_label = QLabel("Selected atoms: (none)")
        layout.addWidget(self.selected_label)
        
        # Fragment preview
        self.fragment_info = QLabel("Fragment 1: (empty)\nFragment 2: (all atoms)")
        layout.addWidget(self.fragment_info)
        
        # Control buttons
        button_layout = QHBoxLayout()
        
        self.read_selection_button = QPushButton("Read Selection")
        self.read_selection_button.clicked.connect(self.read_jmol_selection)
        button_layout.addWidget(self.read_selection_button)
        
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        self.reset_button = QPushButton("Reset View")
        self.reset_button.clicked.connect(self.reset_jmol_view)
        button_layout.addWidget(self.reset_button)
        
        layout.addLayout(button_layout)
        
        # Dialog buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        self.setLayout(layout)
        
    def setup_jmol(self) -> None:
        """Setup Jmol for native selection picking."""
        if not self.jmol:
            self.selected_label.setText("Error: No Jmol connection available")
            return
            
        try:
            # Save current orientation
            if hasattr(self.jmol, 'store_orientation'):
                self.jmol.store_orientation()
            self.saved_orientation = getattr(self.jmol, 'orient', None)
            
            # Save current selectionHalos setting and enable it for visual feedback
            response = self.jmol.send_cmd("show selectionHalos")
            if response:
                for reply_dict in response:
                    reply = reply_dict.get('reply', '')
                    if reply.startswith('ECHO:selectionHalos'):
                        # Extract the current setting (ON or OFF) from "ECHO:selectionHalos ON/OFF"
                        if 'ON' in reply:
                            self.saved_selection_halos = 'ON'
                        else:
                            self.saved_selection_halos = 'OFF'
                        break
                else:
                    self.saved_selection_halos = 'OFF'  # Default assumption
            else:
                self.saved_selection_halos = 'OFF'  # Default assumption
            
            print(f"Saved selectionHalos setting: {self.saved_selection_halos}")
            
            # Clear scene but keep molecule, and clear any selection
            self.jmol.send_cmd("select all; halo off; color cpk")
            self.jmol.send_cmd("select none")  # Start with no atoms selected
            
            # Enable selectionHalos for visual feedback during selection
            self.jmol.send_cmd("selectionHalos ON")
            
            # Setup native Jmol selection mode
            self.jmol.send_cmd("set picking ATOM")
            # Let Jmol handle selection natively - no callbacks needed
            self.selected_label.setText("Ready: Click atoms to select, then press 'Read Selection'")
            
        except Exception as e:
            print(f"Error setting up Jmol for fragment definition: {e}")
            self.selected_label.setText(f"Error: {e}")
    
    def read_jmol_selection(self) -> None:
        """Read the current selection from Jmol and update fragments."""
        if not self.jmol:
            print("Error: No Jmol connection available")
            return
            
        try:
            # First get the number of selected atoms
            response = self.jmol.send_cmd("print {selected}.size")
            if not response:
                print("No response from Jmol")
                return
            
            # Extract the count
            count = 0
            for reply_dict in response:
                reply = reply_dict.get('reply', '')
                if reply.startswith('ECHO:'):
                    count = int(reply[5:])
                    break
            
            print(f"Number of selected atoms: {count}")
            
            if count == 0:
                print("No atoms currently selected in Jmol")
                self.selected_label.setText("No atoms selected - click atoms in Jmol first, then press this button")
                return
            
            # Get each selected atom by index
            selected_atoms = []
            for i in range(count):
                atom_response = self.jmol.send_cmd(f"print {{selected}}[{i}].atomno")
                if atom_response:
                    for reply_dict in atom_response:
                        reply = reply_dict.get('reply', '')
                        if reply.startswith('ECHO:'):
                            atom_num = int(reply[5:]) - 1  # Convert to 0-based indexing
                            selected_atoms.append(atom_num)
                            break
            
            print(f"Parsed selected atoms (0-based): {selected_atoms}")
            
            if selected_atoms:
                # Remove duplicates and sort
                self.selected_atoms = sorted(list(set(selected_atoms)))
                
                # Update visual feedback in Jmol - add halos to selected atoms
                atom_list = ','.join(str(i + 1) for i in self.selected_atoms)
                self.jmol.send_cmd(f"select atomno={atom_list}; halo on; color red")
                self.jmol.send_cmd("select none")  # Clear selection for further picking
                
                self.update_display()
                print(f"Successfully read {len(self.selected_atoms)} atoms: {[x+1 for x in self.selected_atoms]}")
            else:
                print("Could not find any selected atoms")
                self.selected_label.setText("Could not read selection - try selecting atoms again")
                
        except Exception as e:
            print(f"Error reading Jmol selection: {e}")
            import traceback
            traceback.print_exc()
            self.selected_label.setText(f"Error reading selection: {e}")

    def update_display(self) -> None:
        """Update the display with current selection."""
        # Update selected atoms label
        if self.selected_atoms:
            atom_nums = [i + 1 for i in sorted(self.selected_atoms)]  # Convert to 1-based
            self.selected_label.setText(f"Selected atoms: {atom_nums}")
        else:
            self.selected_label.setText("Selected atoms: (none)")
        
        # Update fragment info
        if self.selected_atoms:
            remaining = list(range(self.natoms))
            for atom in self.selected_atoms:
                if atom in remaining:
                    remaining.remove(atom)
            
            frag1_nums = [i + 1 for i in sorted(self.selected_atoms)]
            frag2_nums = [i + 1 for i in sorted(remaining)]
            
            self.fragment_info.setText(
                f"Fragment 1: {frag1_nums} ({len(frag1_nums)} atoms)\n"
                f"Fragment 2: {frag2_nums} ({len(frag2_nums)} atoms)"
            )
        else:
            self.fragment_info.setText(
                f"Fragment 1: (empty)\n"
                f"Fragment 2: all atoms (1-{self.natoms})"
            )

    def clear_selection(self) -> None:
        """Clear all selected atoms."""
        if self.jmol:
            try:
                # Remove all halos
                self.jmol.send_cmd("select all; halo off")
                # Reset to CPK colors
                self.jmol.send_cmd("color cpk")
            except Exception as e:
                print(f"Error clearing selection: {e}")
        
        self.selected_atoms.clear()
        self.update_display()
    
    def reset_jmol_view(self) -> None:
        """Reset Jmol view to original orientation."""
        if self.jmol and self.saved_orientation:
            try:
                self.jmol.send_cmd(self.saved_orientation)
            except Exception as e:
                print(f"Error resetting view: {e}")
    
    def get_fragments(self) -> List[List[int]]:
        """Get the defined fragments (0-based indices).
        
        Returns:
            List of fragment definitions, each as a list of atom indices
        """
        if not self.selected_atoms:
            # If no atoms selected, return empty list - don't assume all atoms in fragment
            print("DEBUG: No atoms selected, returning empty fragment list")
            return []
        
        fragments = []
        
        # Fragment 1: selected atoms
        fragments.append(sorted(self.selected_atoms))
        
        # Fragment 2: remaining atoms
        remaining = []
        for i in range(self.natoms):
            if i not in self.selected_atoms:
                remaining.append(i)
        
        if remaining:
            fragments.append(remaining)
        
        print(f"DEBUG: Returning fragments: {fragments}")
        return fragments
    
    def cleanup_jmol(self) -> None:
        """Clean up Jmol state."""
        # Stop event timer
        if hasattr(self, 'event_timer') and self.event_timer:
            self.event_timer.stop()
            
        if self.jmol:
            try:
                # Remove halos and restore normal display
                self.jmol.send_cmd("select all; halo off; color cpk")
                # Clear any current selection
                self.jmol.send_cmd("select none")
                
                # Restore the original selectionHalos setting
                if self.saved_selection_halos is not None:
                    self.jmol.send_cmd(f"selectionHalos {self.saved_selection_halos}")
                    print(f"Restored selectionHalos setting to: {self.saved_selection_halos}")
                
            except Exception as e:
                print(f"Error cleaning up Jmol: {e}")
    
    def accept(self) -> None:
        """Accept the dialog and clean up."""
        self.cleanup_jmol()
        super().accept()
    
    def reject(self) -> None:
        """Reject the dialog and clean up."""
        self.cleanup_jmol()
        super().reject()
    
    def closeEvent(self, event) -> None:
        """Handle dialog close event."""
        self.cleanup_jmol()
        super().closeEvent(event)


class ColorButton(QtWidgets.QPushButton):
    '''
    Custom Qt Widget to show a chosen color.

    Left-clicking the button shows the color-chooser, while
    right-clicking resets the color to None (no-color).
    '''

    colorChanged = QtCore.Signal(object)

    def __init__(self, *args, color=None, **kwargs):
        super(ColorButton, self).__init__(*args, **kwargs)

        self._color = None
        self._default = color
        self.pressed.connect(self.onColorPicker)

        # Set the initial/default state.
        self.setColor(self._default)

    def setColor(self, color):
        if color != self._color:
            self._color = color
            self.colorChanged.emit(color)

        if self._color:
            self.setStyleSheet("background-color: %s;" % self._color)
        else:
            self.setStyleSheet("")

    def color(self):
        return self._color

    def onColorPicker(self):
        '''
        Show color-picker dialog to select color.

        Qt will use the native dialog by default.

        '''
        dlg = QtWidgets.QColorDialog(self)
        dlg.setStyleSheet("background-color : lightgrey;")
        if self._color:
            dlg.setCurrentColor(QtGui.QColor(self._color))

        if dlg.exec():
            self.setColor(dlg.currentColor().name())

    def mousePressEvent(self, e):
        if e.button() == QtCore.Qt.RightButton:
            self.setColor(self._default)

        return super(ColorButton, self).mousePressEvent(e)

class SaveSpectraDialog(QDialog):
    """Dialog for saving spectral data to a file."""
    
    def __init__(self, parent: Optional[QWidget] = None, 
                 default_path: str = str(WORKDIR), 
                 default_filename: str = "spectrum.txt") -> None:
        super().__init__(parent)
        
        self.setWindowTitle("Save Spectral Data")
        self.setMinimumSize(400, 150)
        
        layout = QVBoxLayout()
        
        # File path input
        path_layout = QHBoxLayout()
        path_label = QLabel("File Path:", self)
        self.path_input = QLineEdit(self)
        self.path_input.setText(default_path)
        browse_button = QPushButton("Browse", self)
        browse_button.clicked.connect(self.browse_file)
        
        path_layout.addWidget(path_label)
        path_layout.addWidget(self.path_input)
        path_layout.addWidget(browse_button)
        
        layout.addLayout(path_layout)
        
        # Filename input
        filename_layout = QHBoxLayout()
        filename_label = QLabel("Filename:", self)
        self.filename_input = QLineEdit(self)
        self.filename_input.setText(default_filename)
        
        filename_layout.addWidget(filename_label)
        filename_layout.addWidget(self.filename_input)
        
        layout.addLayout(filename_layout)
        
        # Dialog buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        self.setLayout(layout)
    
    def browse_file(self) -> None:
        """Open file dialog to select file path."""
        file_dialog = QtWidgets.QFileDialog(self)
        file_dialog.setFileMode(QtWidgets.QFileDialog.Directory)
        if file_dialog.exec():
            selected_dirs = file_dialog.selectedFiles()
            if selected_dirs:
                self.path_input.setText(selected_dirs[0])
    
    def get_file_info(self) -> Dict[str, str]:
        """Get the selected file path and name."""
        return {
            'path': self.path_input.text(),
            'filename': self.filename_input.text()
        }
