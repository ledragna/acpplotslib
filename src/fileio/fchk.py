import os
from math import ceil
import numpy as np
from calc.core_functionality import molecule

NCOLS_FCHK_R = 5
NCOLS_FCHK_I = 6
B2A = 0.52917720859
CM2AU = 1./219474.631371

def fchk_parser(fname):
    """Function that parse the Gaussian Formated checkpoint file and get the
    data required by acp analysis

    Arguments:
        fname {str} -- the file name of the .fchk file
    """
    keys = {'nat': 'Number of atoms',
            'crd': 'Current cartesian coordinates',
            'ian': 'Atomic numbers',
            'atm': 'Vib-AtMass',
            'evc': 'Vib-Modes',
            'apt': 'Dipole Derivatives',
            'aat': 'AAT',
            've2': 'Vib-E2'
            }
    data = molecule()
    name = os.path.split(fname)[-1][:-5]
    data.name = name
    qtt = len(keys)
    with open(fname, 'r') as fopen:
        deriv = 'num derivs'
        line = fopen.readline()
        while line:
            if line.startswith(keys['nat']):
                data.natoms = int(line.split()[-1])
                qtt -= 1
            elif line.startswith(keys['crd']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x)*B2A for x in line.split()])
                tmp = np.array(tmp)
                data.crd = tmp.reshape(-1, 3)
                qtt -= 1
            elif line.startswith(keys['ian']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_I)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([int(x) for x in line.split()])
                data.atnum = np.array(tmp)
                qtt -= 1
            elif line.startswith(keys['atm']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x) for x in line.split()])
                data.atmas = np.array(tmp)
                qtt -= 1
            elif line.startswith(keys['apt']) and deriv not in line:
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x) for x in line.split()])
                tmp = np.array(tmp)
                data.apt = tmp.reshape(-1, 3, 3)
                qtt -= 1
            elif line.startswith(keys['aat']) and deriv not in line:
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(e) for e in line.split()])
                tmp = np.array(tmp)
                data.aat = tmp.reshape(-1, 3, 3)
                qtt -= 1
            elif line.startswith(keys['evc']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(x) for x in line.split()])
                tmp = np.array(tmp)
                data.evec = tmp.reshape(-1, 3*data.natoms)
                qtt -= 1
            elif line.startswith(keys['ve2']):
                nval = int(line.split()[-1])
                nline = ceil(nval/NCOLS_FCHK_R)
                tmp = []
                for _ in range(nline):
                    line = fopen.readline()
                    tmp.extend([float(e) for e in line.split()])
                data.frq = np.array(tmp[:data.nmnum])
                data.rmas = np.array(tmp[data.nmnum:2*data.nmnum])
                qtt -= 1
            if qtt == 0:
                break
            line = fopen.readline()

    return data