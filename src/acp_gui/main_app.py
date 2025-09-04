#!/usr/bin/env python3
"""Main ACP Application - contains the application logic and main window."""

import sys
import os
from typing import Optional, Dict, Any, List
import numpy as np

from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
    QGridLayout,
    QSpacerItem,
    QSizePolicy,
    QFileDialog,
    QStatusBar,
    QDialog,
)
from PySide6.QtGui import QIcon, QAction

# Import our components
from acp_gui.gui_components import (
    WidgetPlot,
    FragmentDialog,
    PropertySelector,
    SpectralParametersWidget,
    ControlButtonsWidget,
    InteractiveFragmentDialog,
    SaveSpectraDialog,
)
from utils.helpers import random_colors

# Import dependencies with error handling
try:
    from calc.core_functionality import ELEMENTS
    from graph.plots import getcolorshexa
    from jmol_interface.wrapper import JmolWrapper
    from fileio.json import read_json
    from fileio.fchk import fchk_parser
    from utils.helpers import readlistnum

    DEPENDENCIES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Some core dependencies not available: {e}")
    DEPENDENCIES_AVAILABLE = False


class ACPMainWindow(QMainWindow):
    """Main window for the ACP application."""

    def __init__(self) -> None:
        QMainWindow.__init__(self)

        # Window properties
        self.title = "ACP Plots"
        self.left = 10
        self.top = 10
        self.width = 1920
        self.height = 1080

        # Data properties
        self.molsdata: Optional[Any] = None
        self.fragments: Optional[List[List[int]]] = None
        self.fname: Optional[str] = None
        self.jmol: Optional[Any] = None
        self.__tmpdata: Optional[Any] = None
        self.__color: List[str] = []
        self.__crstatus: Dict[str, Any] = {"selnm": None, "start": None, "end": None}

        self.setup_window()
        self.setup_menu()
        self.setup_ui()

        # Initialize Jmol if dependencies are available
        if DEPENDENCIES_AVAILABLE:
            try:
                self.jmol = JmolWrapper()
            except Exception as e:
                print(f"Warning: Could not initialize Jmol: {e}")
                self.jmol = None

    def setup_window(self) -> None:
        """Setup the main window properties."""
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Create status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

    def setup_menu(self) -> None:
        """Setup the main menu bar."""
        main_menu = self.menuBar()
        main_menu.setNativeMenuBar(False)

        # File menu
        file_menu = main_menu.addMenu("File")
        # save_spec = main_menu.addMenu('Save Spectra', inactive=True)
        help_menu = main_menu.addMenu("Help")

        # Open action
        open_button = QAction(QIcon(""), "Open", self)
        open_button.setShortcut("Ctrl+O")
        open_button.setStatusTip("Open a file")
        open_button.triggered.connect(self.open_file)
        file_menu.addAction(open_button)

        # Save Spectra action
        saveButton = QAction(QIcon(""), "Save Spectra", self)
        saveButton.setShortcut("Ctrl+S")
        saveButton.setStatusTip("Save a dat file with convoluted spectra")
        saveButton.triggered.connect(self._save_spectra)
        file_menu.addAction(saveButton)

        # Exit action
        exit_button = QAction(QIcon("exit24.png"), "Exit", self)
        exit_button.setShortcut("Ctrl+Q")
        exit_button.setStatusTip("Exit application")
        exit_button.triggered.connect(self._close_jmol)
        file_menu.addAction(exit_button)

    def setup_ui(self) -> None:
        """Setup the main user interface."""
        widget = QWidget(self)
        self.setCentralWidget(widget)

        main_layout = QVBoxLayout(widget)
        grid = QGridLayout()
        main_layout.addLayout(grid)

        # Setup components
        self.setup_property_selector(grid)
        self.setup_spectral_parameters(main_layout)
        self.setup_control_buttons(grid)
        self.setup_plot_widget(main_layout)

    def setup_property_selector(self, grid: QGridLayout) -> None:
        """Setup the property selector."""
        self.prop = PropertySelector()
        grid.addWidget(self.prop, 0, 1)

    def setup_spectral_parameters(self, layout: QVBoxLayout) -> None:
        """Setup spectral parameters widget."""
        self.spectral_params = SpectralParametersWidget()
        layout.addWidget(self.spectral_params)

    def setup_control_buttons(self, grid: QGridLayout) -> None:
        """Setup control buttons."""
        self.controls = ControlButtonsWidget()

        # Connect button signals
        self.controls.groups_button.clicked.connect(self.open_groups)
        self.controls.define_groups_button.clicked.connect(self.define_groups)
        self.controls.interactive_fragments_button.clicked.connect(
            self.define_fragments_interactive
        )
        self.controls.acp_button.clicked.connect(self.calculate_acp)

        # Add to layout
        control_layout = QHBoxLayout()
        control_layout.addWidget(self.controls.groups_button)
        control_layout.addWidget(self.controls.define_groups_button)
        control_layout.addWidget(self.controls.interactive_fragments_button)
        control_layout.addWidget(self.controls.acp_colors_cb)
        control_layout.addItem(QSpacerItem(100, 10, QSizePolicy.Expanding))

        acp_layout = QHBoxLayout()
        acp_layout.addWidget(self.controls.acp_button)
        acp_layout.addItem(QSpacerItem(100, 10, QSizePolicy.Expanding))

        grid.addLayout(control_layout, 0, 0)
        grid.addLayout(acp_layout, 0, 2)

    def setup_plot_widget(self, layout: QVBoxLayout) -> None:
        """Setup the matplotlib plot widget."""
        self.mpl_plot = WidgetPlot(self)
        layout.addWidget(self.mpl_plot)
        try:
            self.mpl_plot.canvas.mpl_connect("pick_event", self.on_pick)
        except Exception as e:
            print(f"Warning: Could not connect pick event: {e}")

    def open_file(self) -> None:
        """Open a molecular structure file."""
        fname = QFileDialog.getOpenFileName(self, "Select FCHK file", ".", "*.fchk")[0]

        if not fname:
            return

        self.fname = os.path.basename(fname)

        if not DEPENDENCIES_AVAILABLE:
            self.status_bar.showMessage(
                "Dependencies not available - file cannot be processed"
            )
            return

        try:
            moldata = fchk_parser(fname)
            self.controls.define_groups_button.setEnabled(True)

            # Enable properties based on available data
            if hasattr(moldata, "apt") and moldata.apt is not None:
                self.prop.enable_property("IR")
            if hasattr(moldata, "aat") and moldata.aat is not None:
                self.prop.enable_property("VCD")

            if self.fragments is not None:
                self.controls.acp_button.setEnabled(True)

            self.title = f"ACP - {self.fname}"
            self.setWindowTitle(self.title)
            self.molsdata = moldata
            self._write_geom_jmol()

        except Exception as err:
            print(f"Error opening file: {err}")
            self.status_bar.showMessage(f"Error opening file: {err}")

    def _save_spectra(self) -> None:
        """Save the currently displayed spectra to a file."""
        if self.__tmpdata is None:
            self.status_bar.showMessage("No spectral data to save")
            return

        try:
            spec_type = self.prop.currentText().lower() if self.prop.currentText() else "spectrum"
            
            dialog = SaveSpectraDialog(
                self,
                default_path=os.getcwd(),
                default_filename=self.fname.replace(".fchk", f"_{spec_type}_gcp_spectra.dat")
                if self.fname
                else f"{spec_type}_gcp_spectra.dat",
            )
            if dialog.exec() == QDialog.Accepted:
                file_info = dialog.get_file_info()
                file_path = file_info['path']
                filename = file_info['filename']
                
                if filename and file_path:
                    # Combine path and filename to get the full file path
                    fname = os.path.join(file_path, filename)
                    
                    # Get the spectral parameters to convolute the data
                    params = self.spectral_params.get_parameters()
                    
                    # Convolute the raw spectral data to get the actual plotted spectrum
                    from calc.core_functionality import convolute
                    convoluted = convolute(self.__tmpdata, params["start"], params["end"], params["hwhm"], spec_type)
                    
                    # Extract x and y values from the convoluted spectrum
                    x_vals = convoluted.xvals
                    y_vals = convoluted.yvals
                    
                    with open(fname, "w") as fopen:
                        header = "# {:10s}".format("Frequency")
                        header += "{:15s}".format("Total")
                        header += "".join(
                            "{:15s}".format(f"Frag{i+1}") for i in range(len(y_vals) - 1)
                        )
                        fopen.write(header + "\n")

                        for i in range(len(x_vals)):
                            line = f"  {x_vals[i]:10.3f} "
                            for y in y_vals:
                                line += f"{y[i]:15.6E}"
                            fopen.write(line + "\n")
                    self.status_bar.showMessage(f"Spectra saved to {fname}")
                else:
                    self.status_bar.showMessage("No filename or path provided")
        except Exception as e:
            print(f"Error saving spectra: {e}")
            self.status_bar.showMessage(f"Error saving spectra: {e}")

    def open_groups(self) -> None:
        """Open a JSON file containing fragment definitions."""
        jname = QFileDialog.getOpenFileName(self, "Select JSON file", ".", "*.json")[0]

        if not jname or not DEPENDENCIES_AVAILABLE:
            return

        try:
            jdata = read_json(jname)
            res = []

            if self.molsdata is not None:
                natm = self.molsdata.natoms
            else:
                natm = int(jdata["molecule"]["natoms"])

            for i in jdata["molecule"]["frags"]:
                if isinstance(i, str):
                    tmp = readlistnum(i["fr_index"], natm)
                else:
                    tmp = i["fr_index"]
                res.append([x - 1 for x in tmp])

            self.fragments = res
            self.__color = random_colors(len(res))

            if self.fname is not None:
                self.controls.acp_button.setEnabled(True)

        except Exception as err:
            print(f"Error opening groups: {err}")
            self.status_bar.showMessage(f"Error opening groups: {err}")

        if self.jmol:
            try:
                self.jmol.send_cmd(self._highlight_fragments())
            except Exception as e:
                print(f"Warning: Could not highlight fragments in Jmol: {e}")

    def define_groups(self) -> None:
        """Define molecular groups through a dialog."""
        if not self.fragments:
            self.fragments = [[0, 1], [2, 3]]  # Default fragments
            self.__color = random_colors(2)

        try:
            colors = FragmentDialog(nfrag=len(self.fragments), colors=self.__color)
            colors.exec()

            if colors.result():
                self.__color = colors._colors

                if self.jmol:
                    try:
                        self.jmol.send_cmd(self._highlight_fragments())
                    except Exception as e:
                        print(f"Warning: Could not highlight fragments: {e}")

                if self.__tmpdata is not None:
                    params = self.spectral_params.get_parameters()
                    spectype = self.prop.currentText().lower()
                    self.mpl_plot.plot(
                        self.__tmpdata,
                        params["start"],
                        params["end"],
                        params["hwhm"],
                        spectype,
                        colors=self.__color,
                    )

        except Exception as e:
            print(f"Error defining groups: {e}")

    def define_fragments_interactive(self) -> None:
        """Define molecular fragments interactively using Jmol atom picking."""
        if not self.jmol:
            print("Error: No Jmol connection available")
            return

        if not self.molsdata:
            print("Error: No molecular data available")
            return

        try:
            natoms = self.molsdata.natoms

            dialog = InteractiveFragmentDialog(
                parent=self,
                jmol=self.jmol,
                natoms=natoms,
                existing_fragments=self.fragments,
            )

            if dialog.exec() == QDialog.Accepted:
                new_fragments = dialog.get_fragments()
                print(f"DEBUG: Dialog returned fragments: {new_fragments}")

                if new_fragments:
                    self.fragments = new_fragments
                    self.__color = random_colors(len(self.fragments))

                    # Enable ACP calculation button now that fragments are defined
                    self.controls.acp_button.setEnabled(True)

                    # Update the plot if we have data - always replot to show new colors
                    if self.__tmpdata is not None:
                        self.calculate_acp()

                    # Highlight fragments in Jmol - this needs to happen after dialog cleanup
                    if self.jmol:
                        try:
                            highlight_cmd = self._highlight_fragments()
                            print(f"DEBUG: Sending highlight command: {highlight_cmd}")
                            self.jmol.send_cmd(highlight_cmd)
                        except Exception as e:
                            print(f"Warning: Could not highlight fragments: {e}")

                    print(f"Fragments defined: {self.fragments}")
                else:
                    print("No fragments were defined - user didn't select any atoms")

        except Exception as e:
            print(f"Error in interactive fragment definition: {e}")

    def calculate_acp(self) -> None:
        """Calculate and display ACP spectra."""
        if not DEPENDENCIES_AVAILABLE or not self.molsdata:
            self.status_bar.showMessage(
                "Cannot calculate - missing dependencies or data"
            )
            return

        params = self.spectral_params.get_parameters()
        prop_cur = self.prop.currentText().lower()

        if not prop_cur:
            self.status_bar.showMessage("Please select a property first")
            return

        try:
            self.molsdata.frag = self.fragments if self.fragments else [[0], [1]]

            # Check if parameters changed
            if (
                params["start"] == self.__crstatus["start"]
                and params["end"] == self.__crstatus["end"]
            ):
                if self.__crstatus["selnm"] and self.controls.acp_colors_cb.isChecked():
                    self._color_acp(self.__crstatus["selnm"], prop_cur)
            else:
                if self.jmol:
                    try:
                        self.jmol.send_cmd(
                            "select all; vector off; vibration off; color cpk"
                        )
                    except Exception as e:
                        print(f"Warning: Jmol command failed: {e}")
                self.__crstatus.update(params)
                self.__crstatus["selnm"] = None

            datas = self.molsdata.frag_spect(prop_cur)
            self.__tmpdata = datas
            self.mpl_plot.plot(
                datas,
                params["start"],
                params["end"],
                params["hwhm"],
                prop_cur,
                colors=self.__color if self.__color else ["blue", "red"],
            )

        except Exception as e:
            print(f"Error calculating ACP: {e}")
            self.status_bar.showMessage(f"Error calculating ACP: {e}")

    def on_pick(self, event) -> None:
        """Handle pick events on the plot."""
        if not DEPENDENCIES_AVAILABLE:
            return

        try:
            datadict = {"vcd": ["rs", "10-44 esu2cm2"], "ir": ["ds", "10-40 esu2cm2"]}

            this = event.ind
            self.__crstatus["selnm"] = this[0]
            prop_cur = self.prop.currentText().lower()

            if not prop_cur or prop_cur not in datadict:
                return

            prop = self.molsdata.get_strg(datadict[prop_cur][0])[this][0]

            status_msg = f"Normal mode: {this[0] + 1:4d} {self.prop.currentText()} {prop:8.4f} {datadict[prop_cur][1]}"
            self.status_bar.showMessage(status_msg)

            self._display_vibration(this[0])

            if self.controls.acp_colors_cb.isChecked():
                self._color_acp(this[0], prop_cur)
        except Exception as e:
            print(f"Error handling pick event: {e}")

    def _write_geom_jmol(self) -> str:
        """Write molecular geometry to Jmol."""
        if not self.jmol or not self.molsdata:
            return ""

        try:
            app = "model"
            modelname = f"{self.fname}"
            str2jmol = f"load data '{app} {modelname}'|"
            str2jmol += f"{self.molsdata.natoms}|"
            str2jmol += f"{modelname}|"

            tmpl = "{a:4s}{b[0]:12.5f}{b[1]:12.5f}{b[2]:12.5f}|"
            for i in range(self.molsdata.natoms):
                str2jmol += tmpl.format(
                    a=ELEMENTS[self.molsdata.atnum[i]], b=self.molsdata.crd[i, :]
                )

            str2jmol += (
                f"end '{app} {modelname}';show data; set background 'lightgrey';"
            )
            self.jmol.send_cmd(str2jmol)
            return str2jmol
        except Exception as e:
            print(f"Error writing geometry to Jmol: {e}")
            return ""

    def _highlight_fragments(self) -> str:
        """Generate Jmol command to highlight fragments."""
        if not self.fragments or not self.__color:
            return ""

        cmd_jmol = ""
        for i, frag in enumerate(self.fragments):
            tmp = "select ({{{}}}); halos on; color halos '{}';"
            atoms = " ".join(str(j) for j in frag)
            color = self.__color[i] if i < len(self.__color) else "#FF0000"
            cmd_jmol += tmp.format(atoms, color)
        return cmd_jmol

    def _color_acp(self, ind: int, prop: str) -> None:
        """Color atoms based on ACP values."""
        if not self.jmol or not self.molsdata:
            return

        try:
            wht = self.molsdata.get_strg(f"{prop}acp")[ind, :]
            cols = getcolorshexa(wht)
            ucols = list(set(cols))
            fcols = {x: [] for x in ucols}

            for i, val in enumerate(cols):
                fcols[val].append(i)

            str2jmol = ""
            for key in fcols:
                tmp = "color ({{{}}}) '{}';"
                atoms = " ".join(str(j) for j in fcols[key])
                str2jmol += tmp.format(atoms, key)

            self.jmol.send_cmd(str2jmol)
        except Exception as e:
            print(f"Error coloring ACP: {e}")

    def _display_vibration(self, index: int) -> None:
        """Display vibrational mode in Jmol."""
        if not self.jmol:
            return

        try:
            self.jmol.store_orientation()
            cmd_to_jmol = ""
            cmd_to_jmol += self._write_data_jmol(index)
            cmd_to_jmol += self._highlight_fragments()
            cmd_to_jmol += (
                "select all; vector 0.05; vector scale 4; vibration on; vibration 0.7;"
            )
            cmd_to_jmol += "set background 'lightgrey';"
            if hasattr(self.jmol, "orient") and self.jmol.orient:
                cmd_to_jmol += self.jmol.orient

            self.jmol.send_cmd(cmd_to_jmol)
        except Exception as e:
            print(f"Error displaying vibration: {e}")

    def _write_data_jmol(self, index: int) -> str:
        """Write vibrational data to Jmol."""
        if not self.molsdata:
            return ""

        try:
            app = "model"
            modelname = f"NM: {index + 1:4d} freq:{self.molsdata.frq[index]:8.2f}"
            str2jmol = f"load data '{app} {modelname}'|"
            str2jmol += f"{self.molsdata.natoms}|"
            str2jmol += f"{modelname}|"

            tmpl = "{a:4s}{b[0]:12.5f}{b[1]:12.5f}{b[2]:12.5f} 0 {c[0]:12.5f}{c[1]:12.5f}{c[2]:12.5f}|"
            lxcrt = self.molsdata.evec.reshape(
                self.molsdata.nmnum, self.molsdata.natoms, 3
            )[index, :, :]

            for i in range(self.molsdata.natoms):
                str2jmol += tmpl.format(
                    a=ELEMENTS[self.molsdata.atnum[i]],
                    b=self.molsdata.crd[i, :],
                    c=lxcrt[i, :],
                )

            str2jmol += f"end '{app} {modelname}';show data;"
            return str2jmol
        except Exception as e:
            print(f"Error writing vibration data: {e}")
            return ""

    def _close_jmol(self) -> None:
        """Close Jmol and exit application."""
        if self.jmol:
            try:
                self.jmol.close()
            except Exception as e:
                print(f"Error closing Jmol: {e}")
        self.close()


def create_application() -> QApplication:
    """Create and return the QApplication instance."""
    return QApplication(sys.argv)


def create_main_window() -> ACPMainWindow:
    """Create and return the main window."""
    return ACPMainWindow()
