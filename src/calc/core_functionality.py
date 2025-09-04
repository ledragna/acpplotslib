"""[summary]

Raises:
    ValueError: [description]
    ValueError: [description]
    ValueError: [description]
    ValueError: [description]
    ValueError: [description]

Returns:
    [type] -- [description]
"""
import copy
from math import log, pi
import numpy as np
from typing import List
from utils.gxx_units_np import PhyCon


# Constants for conversion
m2ang = 1.0e10
hbar = PhyCon(4)*m2ang**2/(2*pi*PhyCon(2))
hc = PhyCon(4)*PhyCon(9)*1.0e18
# -- Conversion from mass-weighted to dimensionless normal coords
MWQ2q = 1./(np.sqrt(2*pi*PhyCon(9)/hbar)*PhyCon(1))
# -- Conversion from sqrt(Eh/amu) to cm^-1
eval2cm2 = MWQ2q**2/hc * PhyCon(8)*1.0e18
# -- Electric dipole from au to statC.cm
edip_conv = PhyCon(3)*PhyCon(1)*1.0e-8
# -- Magnetic dipole from au to statA.cm^2
mdip_conv = 1.0e4*PhyCon(4)/(2*pi)*PhyCon(3) / PhyCon(2)
# -- IR to epsilon
IR2EPS = 1.0e-47*8.*pi**3*PhyCon(5) / (3000.*PhyCon(4)*PhyCon(9)*log(10))
# -- VCD to Depsilon
VCD2DE = 1.0e-51*32.*pi**3*PhyCon(5) / (3000.*PhyCon(4)*PhyCon(9)*log(10))
ds_fact = MWQ2q**2 * edip_conv**2*1.0e40/2.
rs_fact = edip_conv*mdip_conv*1.0e44/PhyCon(9)


#typing alias
LLInt = List[List[int]]


NCOLS_FCHK_R = 5
NCOLS_FCHK_I = 6
B2A = 0.52917720859
CM2AU = 1./219474.631371
ELEMENTS = ['',
    'H',                                                                                'He',
    'Li','Be'                                                  ,'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
    'Na','Mg'                                                  ,'Al','Si','P' ,'S' ,'Cl','Ar',
    'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
    'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
        'Hf','Ta',' W','Re','Os','Ir','Pt', 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',
    'Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',
        'Rf','Db','Sg','Bh','Hs','Mt','Ds']






class ccmimate():
    def __init__(self, eng, peaks):
        self.etenergies = eng
        self.obs = peaks


def calc_broad(data, start, end, hwhm):
    """
    Calculates a broadened Gaussian spectrum for the given data within a specified energy range.
    Uses optimized Cython implementation for high-performance spectral broadening.

    Parameters:
        data: An object containing 'etenergies' (energy values) and 'obs' (observable values).
        start (float): The starting energy value for the spectrum.
        end (float): The ending energy value for the spectrum.
        hwhm (float): The half-width at half-maximum (HWHM) for the Gaussian broadening.

    Returns:
        OptimizedSpectrum: An instance representing the broadened spectrum over the specified range.

    Notes:
        - The number of points in the spectrum is set to five times the energy range.
        - Uses Cython-accelerated broadening with OpenMP parallelization.
        - Requires the broadening.pyx module to be compiled.
    """
    from .broadening import broadlor
    import multiprocessing
    
    numpts = int((end - start) * 5)
    
    # Create x-axis for the spectrum
    xaxis = np.linspace(start, end, numpts, dtype=np.double)
    
    # Convert input data to required format
    x_trans = np.array(data.etenergies, dtype=np.double)
    y_trans = np.array(data.obs, dtype=np.double)
    
    # Use optimal number of threads
    num_threads = min(multiprocessing.cpu_count(), 4)
    
    # Call optimized Cython function
    spectrum = broadlor(x_trans, y_trans, xaxis, hwhm, num_threads)
    
    # Create a result object that mimics GaussianSpectrum interface
    class OptimizedSpectrum:
        def __init__(self, xvalues, spectrum_data):
            self.xvalues = xvalues
            self.spectrum = spectrum_data.reshape(1, -1)  # Match GaussianSpectrum format
    
    return OptimizedSpectrum(xaxis, spectrum)

class molecule():
    def __init__(self,
                 name='default',
                 natoms=1,
                 crd=None,
                 atnum=None,
                 atmas=None,
                 mult=1,
                 chrg=0,
                 evec=None,
                 apt=None,
                 aat=None,
                 frq=None,
                 rmas=None):
        self.name = name
        self.natoms = natoms
        self.crd = crd
        self.atnum = atnum
        self.atmas = atmas
        self.mult = mult
        self.chrg = chrg
        self.evec = evec
        self.apt = apt
        self.aat = aat
        self.frq = frq
        self.rmas = rmas
        self.__frag = None
        self.__strg = {}

    @property
    def name(self):
        """
        Return the name of the molecular system
        """
        return self.__name

    @name.setter
    def name(self, val):
        """Sets the name of the molecular system

        Arguments:
            val {str} -- the name of the molecular system
        """
        self.__name = val

    @property
    def natoms(self):
        """
        int: natoms
        """
        return self.__natoms

    @natoms.setter
    def natoms(self, val):
        self.__natoms = val
        try:
            self.__nmnum = val * 3 - 6
        except Exception:
            self.__nmnum = 0

    @property
    def nmnum(self):
        return self.__nmnum

    @property
    def crd(self):
        return self.__crd

    @crd.setter
    def crd(self, val):
        self.__crd = val

    @property
    def atnum(self):
        return self.__atnum

    @atnum.setter
    def atnum(self, val):
        self.__atnum = val

    @property
    def atmas(self):
        return self.__atmas

    @atmas.setter
    def atmas(self, val):
        self.__atmas = val

    @property
    def mult(self):
        return self.__mult

    @mult.setter
    def mult(self, val):
        self.__mult = val

    @property
    def chrg(self):
        return self.__chrg

    @chrg.setter
    def chrg(self, val):
        self.__chrg = val

    @property
    def evec(self):
        return self.__evec

    @evec.setter
    def evec(self, val):
        self.__evec = val

    @property
    def apt(self):
        return self.__apt

    @apt.setter
    def apt(self, val):
        self.__apt = val

    @property
    def aat(self):
        return self.__aat

    @aat.setter
    def aat(self, val):
        self.__aat = val

    @property
    def frq(self):
        return self.__frq

    @frq.setter
    def frq(self, val):
        self.__frq = val

    @property
    def rmas(self):
        return self.__rmas

    @rmas.setter
    def rmas(self, val):
        self.__rmas = val

    @property
    def frag(self) -> LLInt:
        return self.__frag

    @frag.setter
    def frag(self, val: LLInt) -> None:
        tmp_full = [item for sublist in val for item in sublist]
        natms = len(tmp_full)
        unic_atms = set(tmp_full)
        lunic_atms = len(unic_atms)
        if natms != lunic_atms:
            print("WARNING: some atoms present more than one fragments")
        try:
            unic_atms = np.array(list(unic_atms))
            if (unic_atms < 0).any() or (unic_atms >= self.nmnum).any():
                self.__frag = None
                raise ValueError
            else:
                excluded = list(set(range(self.natoms)) - set(tmp_full))
                if excluded:
                    val.append(excluded)
                self.__frag = [tuple(x) for x in val]
        except Exception as err:
            print(err)


    def _calc_lx(self):
        """
        return the cartesian L
        """
        self.__strg['lx'] = self.evec/np.sqrt(self.rmas)[:, np.newaxis]

    def _calc_v(self):
        """
        calc
        vdip = dot(apt, apt.T)
        vrot = dot(apt, aat.T)
        and stores internally
        """

        self.__strg['vdip'] = np.dot(self.apt.reshape(-1, 3),
                                     self.apt.reshape(-1, 3).T)
        self.__strg['vrot'] = np.dot(self.apt.reshape(-1, 3),
                                     self.aat.reshape(-1, 3).T)

    def _calc_p(self):
        """
        calc the invariant
        J = sum_{a,b} L_a * V_{a,b} * L_b
        and stores internally

        """
        if 'lx' not in self.__strg:
            self._calc_lx()
        if 'vdip' not in self.__strg:
            self._calc_v()
        shape = (self.nmnum, self.natoms, 3, self.natoms, 3)
        llt = np.einsum('ki,kj->kij', self.__strg['lx'], self.__strg['lx'])
        tmp_p = (llt * self.__strg['vdip'][np.newaxis, :]).reshape(shape)
        tmp_r = (llt * self.__strg['vrot'][np.newaxis, :]).reshape(shape)
        self.__strg['dipp'] = np.einsum('ijklm->ijl', tmp_p)
        self.__strg['rotp'] = np.einsum('ijklm->ijl', tmp_r)

    def get_strg(self, name):
        return self.__strg[name]

    def _calc_obs(self):
        """
        Calculate dipole strength and rotational strength
        and stores them
        """
        if 'lx' not in self.__strg:
            self._calc_lx()
        if 'vdip' not in self.__strg:
            self._calc_v()
        lxcrt = self.__strg['lx']
        vdip = self.__strg['vdip']
        vrot = self.__strg['vrot']

        self.__strg['ds'] = np.einsum('li,ij,lj->l',
                                      lxcrt, vdip, lxcrt) * ds_fact / self.frq
        self.__strg['rs'] = np.einsum('li,ij,lj->l',
                                      lxcrt, vrot, lxcrt) * rs_fact

    def _calc_gcm(self):
        """
        compute group contribution matrices
        """
        keys = {'ir': ['dip', ds_fact / self.frq],
                'vcd': ['rot', rs_fact * np.ones(self.frq.shape)]
                }
        if self.frag is None:
            print("frag not set")
            raise ValueError
        #l dimension is the normal modes
        if 'dipp' not in self.__strg:
            self._calc_p()
        nfrag = len(self.frag)
        nmodes = len(self.frq)
        natms = self.natoms
        for obs in keys:
            acm = self.get_strg("{}p".format(keys[obs][0]))
            gcm = np.zeros((nmodes, nfrag, nfrag))
            for j in range(nmodes):
                tmp_gcm = np.zeros((nfrag, natms))
                for i in range(nfrag):
                    tmp_gcm[i, :] = acm[j, self.frag[i], :].sum(axis=0)
                for i in range(nfrag):
                    gcm[j, :, i] = tmp_gcm[:, self.frag[i]].sum(axis=1)
                gcm[j, :, :] *= keys[obs][1][j]
            self.__strg['{}gcm'.format(obs)] = gcm

    def _calc_acp(self):
        """
        calc acp and stores it
        """
        keys = {'ir': ['dip', ds_fact / self.frq],
                'vcd': ['rot', rs_fact * np.ones(self.frq.shape)]
                }
        #l dimension is the normal modes
        if 'lx' not in self.__strg:
            self._calc_lx()
        lxcrt = self.__strg['lx'].reshape(self.nmnum, self.natoms, 3)
        if 'dipp' not in self.__strg:
            self._calc_p()
        if 'vdip' not in self.__strg:
            self._calc_v()

        for obs in keys:
            # Not Working
            jobs = self.__strg['{}p'.format(keys[obs][0])]
            mat = self.__strg['v{}'.format(keys[obs][0])].reshape(self.natoms, 3,
                                                               self.natoms, 3)
            diagonal = mat.diagonal(0, 0, 2)
            # Lxa_Vaa = V_a,alpha= \sum_beta Lx_{a,beta} V_{a,beta,a,alpha}
            lxavaa = np.einsum('ijk,lki->ljk', diagonal, lxcrt)
            lxavaa = np.einsum('lij->lji', lxavaa)
            # Vaa_Lxa = V_alpha,a= \sum_beta  V_{a,alpha,a,beta} Lx_{a,beta}
            vaalxa = np.einsum('ijk,lkj->lik', diagonal, lxcrt)
            vaalxa = np.einsum('lij->lji', vaalxa)
            # Lxa_Vab = V_a,b,alpha= \sum_beta Lx_{a,beta} V_{a,beta,b,alpha}
            lxavab = np.einsum('lij,ijkm->likm', lxcrt, mat)
            # Vab_Lxb = V_a,b,alpha= \sum_beta  V_{a,alpha,b,beta} Lx_{b,beta}
            vablxb = np.einsum('lkm,ijkm->lkji', lxcrt, mat)
            # reorder
            vablxb = np.einsum('lijk->lkij', vablxb)
            #Norms
            norm_lxavaa = np.sqrt(np.einsum('lij,lij->li', lxavaa, lxavaa))
            norm_vaalxa = np.sqrt(np.einsum('lij,lij->li', vaalxa, vaalxa))
            del lxavaa, vaalxa
            norm_lxavab = np.sqrt(np.einsum('lijk,lijk->lij', lxavab, lxavab))
            norm_vablxb = np.sqrt(np.einsum('lijk,lijk->lij', vablxb, vablxb))
            del lxavab, vablxb

            #coefficients
            den = norm_lxavab[:, :, :] + norm_lxavaa[:, :, np.newaxis] +\
                  norm_vablxb[:, :, :] + norm_vaalxa[:, np.newaxis, :]
            mask = np.abs(den) < 1e-8
            rab = (norm_lxavab[:, :, :] + norm_lxavaa[:, :, np.newaxis]) / den
            rba = (norm_vablxb[:, :, :] + norm_vaalxa[:, np.newaxis, :]) / den
            rab[mask] = 0.
            rba[mask] = 0.
            self.__strg['{}acp'.format(obs)] = (np.einsum('lij,lij->li', rab, jobs) +\
                                                np.einsum('lij,lij->lj', rba, jobs))
            self.__strg['{}acp'.format(obs)] *= keys[obs][1][:, np.newaxis]

    def frag_spect(self, qtt):
        """
        Calc the spectrum of the fragments
        """
        datadict = {'vcd': ['vcdacp', 'rs', VCD2DE],
                    'ir': ['iracp', 'ds', IR2EPS]
                   }

        result = []
        try:
            acp = self.get_strg(datadict[qtt][0])
        except KeyError:
            self._calc_acp()
            acp = self.get_strg(datadict[qtt][0])
        try:
            obs = self.get_strg(datadict[qtt][1])
        except KeyError:
            self._calc_obs()
            obs = self.get_strg(datadict[qtt][1])
        # the total
        result.append(ccmimate(self.frq, obs * datadict[qtt][2]))
        for frg in self.frag:
            result.append(ccmimate(self.frq, acp[:, frg].sum(axis=1)))

        return result


class conv_spect():
    """Object containing the convoluted spectra on a xvalues grid of points

    Raises:
        ValueError: [description]
        ValueError: [description]

    Returns:
        [type] -- [description]
    """

    def __init__(self, xvals, yvals):
        self.xvals = xvals
        self.yvals = yvals

    @property
    def xvals(self):
        return self.__xvals

    @xvals.setter
    def xvals(self, val):
        self.__xvals = val

    @property
    def yvals(self):
        return self.__yvals

    @yvals.setter
    def yvals(self, val):
        self.__yvals = val

    def get_ymax(self):
        maxim = []
        for spc in self.yvals:
            maxim.append(np.abs(spc).max())
        return np.array(maxim).max()

    def _same_system(self, other):
        same = False
        if (np.abs(self.xvals.min() - other.xvals.min()) < 1e-6 and
                np.abs(self.xvals.max() - other.xvals.max()) < 1e-6 and
                self.xvals.shape == other.xvals.shape and
                len(self.yvals) == len(other.yvals)):
            same = True
        return same

    def __add__(self, other):
        """
        add two set of spectra if share the same range and have the same
        number of yvalues
        """
        try:
            if not self._same_system(other):
                raise ValueError ('Different spectra')
            tmp = copy.deepcopy(self)
            tmp.yvals = []
            for i, yth in enumerate(self.yvals):
                tmp.yvals.append(yth + other.yvals[i])

            return tmp
        except ValueError as err:
            print(err)

    def __iadd__(self, other):
        """
        overload +=
        """
        try:
            if not self._same_system(other):
                raise ValueError('Different spectra')
            for i, yth in enumerate(other.yvals):
                self.yvals[i] += yth
            return self
        except ValueError as err:
            print(err)

    def __sub__(self, other):
        """
        sub two set of spectra if share the same range and have the same
        number of yvalues
        """
        try:
            if not self._same_system(other):
                raise ValueError('Different spectra')
            tmp = copy.deepcopy(self)
            tmp.yvals = []
            for i, yth in enumerate(self.yvals):
                tmp.yvals.append(yth - other.yvals[i])

            return tmp
        except ValueError as err:
            print(err)

    def __isub__(self, other):
        """
        overload -=
        """
        try:
            if not self._same_system(other):
                raise ValueError('Different spectra')
            for i, yth in enumerate(other.yvals):
                self.yvals[i] -= yth
            return self
        except ValueError as err:
            print(err)

    def __mul__(self, param):
        """
        multiply the data set per parameters
        """
        tmp = copy.deepcopy(self)
        tmp.yvals = []
        for yth in self.yvals:
            tmp.yvals.append(yth * param)

        return tmp

    def __imul__(self, param):
        """
        overload *=
        """
        for yth in self.yvals:
            yth *= param

        return self

    def __truediv__(self, param):
        """
        overload true division
        """
        tmp = copy.deepcopy(self)
        tmp.yvals = []
        for yth in self.yvals:
            tmp.yvals.append(yth / param)

        return tmp


def convolute(data, start, end, hwhm, spectype):
    """
    Return an object with x vals and spect
    """
    datadict = {'vcd': VCD2DE,
                'ir': IR2EPS
                   } 
    if spectype not in ['ir', 'vcd']:
        raise ValueError('spectype not recognized')

    smin = start
    smax = end
    # BUG to check
    hwhmo = hwhm
    numpts = int((smax-smin)*5)
    xvalues = np.arange(numpts)*float(end - start)/(numpts-1) + start
    result = []
    for i, subdata in enumerate(data):
        # Non so se tenere un hwhm diverso per i frammenti
        hwhm_n = hwhmo # if not i else 4
        values = calc_broad(subdata, smin, smax, hwhm_n)
        result.append(values.spectrum[0, :]*datadict[spectype]*xvalues)
    return conv_spect(xvalues, result)


def frag_spect(mol, frag, qtt):
    """
    Calc the spectrum of the fragments
    """
    datadict = {'vcd': ['vcdacp', 'rs'],
                'ir': ['iracp', 'ds']
               }
    result = []
    acp = mol.get_strg(datadict[qtt][0])
    obs = mol.get_strg(datadict[qtt][1])
    # the total
    result.append(ccmimate(mol.frq, obs))
    for frg in frag:
        result.append(ccmimate(mol.frq, acp[:, frg].sum(axis=1)))
    return result

def calc_gcm(mol, frag=None, qtt='vcd'):
    """
    Calc the spectrum of the fragments
    """
    datadict = {'vcd': ['rotp'],
                'ir': ['dipp']
               }
    if frag is None:
        frag = mol.frag
    nfrag = len(frag)
    acm = mol.get_strg(datadict[qtt][0])
    nmodes = acm.shape[0]
    natms = acm.shape[1]
    gcm = np.zeros((nmodes, nfrag, nfrag))
    for j in range(nmodes):
        tmp_gcm = np.zeros((nfrag, natms))
        for i in range(nfrag):
            tmp_gcm[i, :] = acm[j, frag[i], :].sum(axis=0)
        for i in range(nfrag):
            gcm[j, :, i] = tmp_gcm[:, frag[i]].sum(axis=1)
    return gcm
