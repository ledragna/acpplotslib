from .gxx_exceptions import PhysError

def PhyCon(id):
    ToAng = 0.52917720859
    ToKG = 1.660538782e-27
    ToE = 1.602176487e-19
    Planck = 6.62606896e-34
    Avog = 6.02214179e23
    JPCal = 4.184
    MPerB = 0.0
    Hartre = 4.35974394e-18
    SLight = 2.99792458e10
    Boltz = 1.3806504e-23
    FineSC = 137.035999679
    EMKG = 0.0
    VolMol = 22.413996e-3
    EMM = 928.476377e-26
    PRM = 1.672621637e-27
    GFree = 2.0023193043622

    # Modified values
    ToE    = ToE*SLight/10.
    MPerB  = ToAng/1.0e10
    FineSC = 1./FineSC
    EMKG   = Hartre*1.0e4/(SLight*FineSC)**2

    if id == 0:
        return '''
 1: Angstroms per Bohr
 2: Kg per atomic mass unit
 3: ESU per electron
 4: Planck constant in J.s
 5: Avogadro constant
 6: Joules per calorie
 7: Meters per Bohr
 8: Joules per Hartree
 9: Speed of light in cm/s
10: Boltzman constant in J/K
11: Inverse Fine structure constant
12: Electron mass in Kg == Kg per atomic unit of mass
13: Molar volume of ideal gas in m**3 at 273.15 K
14: Electron Magnetic Moment (J/Tesla) (sign flipped)
15: Proton Rest Mass (Kg)
16: Free-electron g-factor
'''
    elif id == 1:
        return ToAng
    elif id == 2:
        return ToKG
    elif id == 3:
        return ToE
    elif id == 4:
        return Planck
    elif id == 5:
        return Avog
    elif id == 6:
        return JPCal
    elif id == 7:
        return MPerB
    elif id == 8:
        return Hartre
    elif id == 9:
        return SLight
    elif id == 10:
        return Boltz
    elif id == 11:
        return FineSC
    elif id == 12:
        return EMKG
    elif id == 13:
        return VolMol
    elif id == 14:
        return EMM
    elif id == 15:
        return PRM
    elif id == 16:
        return GFree
    else:
        raise PhysError(id)
