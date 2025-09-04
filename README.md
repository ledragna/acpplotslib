# acpplots Library

## Overview
The `acpplots` library is designed to visulize group contribution patters in spectroscopic simulations.
It is based on the work of [Hug](https://pyvib2.sourceforge.net/) and which was implemented in the software [PyVib2](https://pyvib2.sourceforge.net/) (written in python2 and no longer maintained).
Based on this partition of the contrinutions, the graphical representation is new and rely on Qt and matplotlib for the spectra and on Jmol for the 3D models.
The `jmol_interface` module provides a simple socket-based interface with Jmol, which can also be used for other purposes.


## Installation
To install the `acpplots` library, clone the repository and run the following command in the project directory:

```
pip install .
```

Alternatively, you can install the required dependencies listed in `requirements.txt`:

```
pip install -r requirements.txt
```

## Usage
After installation, you can import the library and use its functionalities as follows:

```python
from acpplots.acp_gui import gui_components
from acpplots.acp_rev import core_functionality
from acpplots.jmol_interface import wrapper
from acpplots.utils import helpers
```

## TODO
 - [ ] only IR and VCD supported, add support for Raman and ROA
 - [ ] a lot...

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
