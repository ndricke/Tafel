import numpy as np
import scipy as sc



class ElecMech(object):
    """concentrations: dictionary of the concentrations of non-intermediate species (pH, [O2], ect)
    rate_constants: dictionary of the rate constants for every step in the reaction network
    ep: for finding pH/concentration dependence, we must assume an observed constant current with varying concentration
    """
    def __init__(s, k_rate=None, conc=None, dG=0., a=0.5):
        s.f = 38.949
        s.R = 0.008314 #In units of kJ/(mol*K)
        s.T = 298.15
        s.kT = 8.61733*10**-5 * s.T #in eV/K
        s.a = a
        #s.dG = dG
        #s.Keq = np.exp(-dG/(s.R*s.T))
        if k_rate is not None: s.setConsts(k_rate)
        if conc is not None: s.setConc(conc)

    @property
    def pH(self):
        return self._pH

    @pH.setter
    def pH(self, value):
        self.H = 10**(-1.*value)
        self._pH = value

    def PCET(s, ks, V=None, H=None, mech='acid'):
        if mech == 'acid':
            return ks[0]*H*np.exp(-0.5*s.f*V), ks[1]*np.exp(0.5*s.f*V)
        elif mech == 'base':
            OH = (10.**-14)/H
            return ks[0]*np.exp(-0.5*s.f*V), ks[1]*OH*np.exp(0.5*s.f*V)
        elif mech == 'both':
            ka, kna = s.PCET(ks[:2], V, H, mech='acid')
            kb, knb = s.PCET(ks[:2], V, H, mech='base')
            return ka + kb*ks[2], kna + knb*ks[2]

    def irrevPCET(s, k, V=None, H=None):
        return k*H*np.exp(-0.5*s.f*V)

    def PT(s, ks, H=None, mech='acid'):
        return ks[0]*H, ks[1]

    def ET(s, ks, V=None, mech='acid'):
        return ks[0]*np.exp(-0.5*s.f*V), ks[1]*np.exp(0.5*s.f*V)

    def kMarcus(s, dG, lam, A):
        return A*np.exp(-(dG+lam)**2/(4.*lam*s.kT))

    def kMHC(s, ieta, ilam):
        eta, lam = -1.*ieta/(s.kT), ilam/(s.kT)
        sqrt_lam = np.sqrt(lam)
        exp_part = 1/(1 + np.exp(-eta))
        erf_part = 1 - sc.special.erf((lam-np.sqrt(1+sqrt_lam+eta**2))/(2.*sqrt_lam))
        return np.sqrt(np.pi) * sqrt_lam * exp_part * erf_part

    def gen_k_edge(s, dG, TS):
        if dG <= 0.:
            dG_fwd, dG_back = 0, dG
        else:
            dG_fwd, dG_back = dG, 0
        k1, kn1 = np.exp(-(TS+dG_fwd)/s.kT), np.exp(-(TS-dG_back)/s.kT)
        return k1, kn1

    def lograte(s, V):
        J = s.rate(V)
        return np.log10(np.abs(J))

    @staticmethod
    def rev2(k1, k2, kn1, kn2):
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)

    @property
    def pH(self):
        return self._pH

    @pH.setter
    def pH(self, value):
        self.H = 10.**(-value)
        self._pH = value

    ##To be defined in child class
    def rate(s, V):
        pass

    ##To be defined in child class
    def setConsts(s, p):
        pass

    ##To be defined in child class
    def setConc(s, p):
        pass

class ClassicBV(ElecMech):

    def rate(s, V):
        return s.k1*np.exp(-0.5*s.f*V) - s.kn1*np.exp(0.5*s.f*V)

    def setConsts(s, p):
        s.k1, s.kn1 = p

    def genConsts(s, dGi_Ti):
        s.k1, s.kn1 = s.gen_k_edge(dGi_Ti[0], dGi_Ti[1])

    def setConc(s, p):
        s.pH = p

class BV(ElecMech):

    def rate(s, V):
        return s.k1*10**(-s.pH)*np.exp(-0.5*s.f*V) - s.kn1*np.exp(0.5*s.f*V)

    def setConsts(s, p):
        s.k1, s.kn1 = p

    def genConsts(s, dGi_Ti):
        s.k1, s.kn1 = s.gen_k_edge(dGi_Ti[0], dGi_Ti[1])

    def setConc(s, p):
        s.pH = p

##Butler-Volmer, with reversible mechanism
class MultiBV(ElecMech):
    def rate(s, V):
        fwd = (s.A*10.**(-s.pH) + s.B)*np.exp(-s.a*V*s.f)
        bck = (s.A + s.B*10.**(s.pH))*np.exp((1.-s.a)*V*s.f)*s.C
        return fwd - bck

    def setConsts(s, p):
        s.A, s.B, s.C = p

    def setConc(s, p):
        s.pH = p

class MHC(ElecMech):

    def rate(s, V):
        eta = V + s.dG
        k1, kn1 = s.kMHC(eta, s.lam), s.kMHC(-eta, s.lam)
        return s.A*(k1 - kn1)

    def setConsts(s, p):
        s.A, s.dG, s.lam  = p

    def setConc(s, p):
        s.pH = p

"""Mechanism with a PCET followed by split proton and electron transfer"""
class Rev2PTET(ElecMech):

    def rate(s, V):
        H = 10.**(-s.pH)
        k1, kn1 = s.PCET((s.k1,s.kn1), V=V, H=H, mech='acid')
        k2, kn2 = s.PT((s.k2,s.kn2), H=H, mech='acid')
        k3, kn3 = s.ET((s.k3,s.kn3), V=V)
        return (k1*k2*k3 - kn1*kn2*kn3)/(k1*k2 + k1*k3 + k2*k3 + k3*kn1 + k1*kn2 + kn1*kn2 + k2*kn3 + kn1*kn3 + kn2*kn3)

    def setConsts(s, p):
        s.k1, s.k2, s.k3, s.kn1, s.kn2, s.kn3 = p

    def setConc(s, p):
        s.pH = p

"""Cycle with 2 reversible eletron transfers, one as a PCET, one as O2-CET"""
class Rev2PcetO2(ElecMech):

    def rate(s, V):
        H = 10.**(-s.pH)
        k1, kn1 = s.PCET((s.k1,s.kn1), V=V, H=H, mech='acid')
        k2, kn2 = s.PT((s.k2,s.kn2), H=s.O2, mech='acid')
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)

    def setConsts(s, p):
        s.k1, s.k2, s.kn1, s.kn2 = p

    def setConc(s, p):
        s.pH = p[0]
        s.O2 = p[1]


"""Cycle with 2 reversible, acid-mediated PCET steps"""
class RevPcet2(ElecMech):

    def rate(s, V):
        H = 10.**(-s.pH)
        k1, k2 = s.k1*H*np.exp(-0.5*s.f*V), s.k2*H*np.exp(-0.5*s.f*V)
        kn1, kn2 = s.kn1*np.exp(0.5*s.f*V), s.kn2*np.exp(0.5*s.f*V)
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)

    def setConsts(s, p):
        s.k1, s.k2, s.kn1, s.kn2 = p

    def setConc(s, p):
        s.pH = p

    """Generate the rate constants for this reaction from intermediate and TS energies"""
    def genConsts(s, dGi_Ti):
        s.k1, s.kn1 = s.gen_k_edge(dGi_Ti[0], dGi_Ti[2])
        s.k2, s.kn2 = s.gen_k_edge(dGi_Ti[1], dGi_Ti[3])

class RevPcet2MHC(ElecMech):

    def rate(s, V):
        eta1, eta2 = V + s.dG1, V + s.dG2
        k1, kn1 = s.kMHC(eta1, s.lam1), s.kMHC(-eta1, s.lam1)
        k2, kn2 = s.kMHC(eta2, s.lam2), s.kMHC(-eta2, s.lam2)
        return s.A * s.rev2(k1, k2, kn1, kn2)

    def setConsts(s, p):
        s.A, s.dG1, s.dG2, s.lam1, s.lam2 = p

    def setConc(s, p):
        s.pH = p


class RevPcet2ab(ElecMech):

    def rate(s, V):
        H = 10.**(-s.pH)
        k1, kn1 = s.PCET((s.k1,s.kn1, s.kab1), V=V, H=H, mech='both')
        k2, kn2 = s.PCET((s.k2,s.kn2, s.kab2), V=V, H=H, mech='both')
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)

    def setConsts(s, p):
        s.k1, s.k2, s.kn1, s.kn2, s.kab1, s.kab2 = p

    def setConc(s, p):
        s.pH = p

    def genConsts(s, dGi_Ti):
        """
        k1, kn1: PCET step with both acidic and basic mechanisms
        k2, kn2: PCET step with both acidic and basic mechanisms
        kab1: ratio of basic/acidic mechanism barrier for step 1
        kab2: ratio of basic/acidic mechanism barrier for step 2
        """
        s.k1, s.kn1 = s.gen_k_edge(dGi_Ti[0], dGi_Ti[2])
        s.k2, s.kn2 = s.gen_k_edge(dGi_Ti[1], dGi_Ti[3])
        s.kab1 = np.exp(-dGi_Ti[4] + dGi_Ti[2])
        s.kab2 = np.exp(-dGi_Ti[5] + dGi_Ti[3])


"""Cycle with 2 reversible, acid- and base-mediated PCET steps"""
class Rev2PTETab(ElecMech):

    def rate(s, V):
        H = 10.**(-s.pH)
        k1, kn1 = s.PCET((s.k1,s.kn1, s.kab1), V=V, H=H, mech='both')
        k2, kn2 = s.PT((s.k2,s.kn2, s.kab2), H=H, mech='both')
        k3, kn3 = s.ET((s.k3,s.kn3), V=V)
        return (k1*k2*k3 - kn1*kn2*kn3)/(k1*k2 + k1*k3 + k2*k3 + k3*kn1 + k1*kn2 + kn1*kn2 + k2*kn3 + kn1*kn3 + kn2*kn3)

    def setConsts(s, p):
        s.k1, s.k2, s.k3, s.kn1, s.kn2, s.kn3, s.kab1, s.kab2 = p

    def setConc(s, p):
        s.pH = p

"""Cycle with a reversible, pH dependent trap"""
class Cyc3Trap(ElecMech):

    def rate(s, V):
        H = 10.**(-s.pH)
        k1, kn1 = s.PCET((s.k1,s.kn1), V=V, H=H, mech='acid')
        k2 = s.irrevPCET(s.k2, V=V, H=H)
        k3 = s.irrevPCET(s.k3, V=V, H=H)
        ktr, kntr = s.PCET((s.ktr,s.kntr), V=V, H=H, mech='acid')
        return k1*k2*k3*kntr/(k1*k3*ktr + k1*k2*kntr + k1*k3*kntr + k2*k3*kntr + k3*kn1*kntr)

    def setConsts(s, p):
        s.k1, s.kn1, s.k2, s.k3, s.ktr, s.kntr = p

    def genConsts(s, dGi_Ti):
        s.k1, s.kn1 = s.gen_k_edge(dGi_Ti[0], dGi_Ti[1])
        s.k2 = np.exp(-dG_Ti[2]/s.kT)
        s.k3 = np.exp(-dG_Ti[3]/s.kT)
        s.ktr, s.kntr = s.gen_k_edge(dGi_Ti[4], dGi_Ti[5])

    def setConc(s, p):
        s.pH = p




if __name__ == "__main__":
    import SimTafel
    import matplotlib.pyplot as plt

    dom = 400
    pH_dom = np.linspace(1,14,dom)
    V_dom = np.linspace(-1,1,dom)

    mech = BV(conc=7.) #Can create regime's with 120mV/dec transition to 40mV/dec
    mech.genConsts((-0.5,0.002))

    #mech = RevPcet2(conc=3.)
    #mech.genConsts((-0.5,0.2, 0.3, 0.1))

    sim = SimTafel.SimTafel(mech)

    #sim.plotTafel(V_dom)
    sim.plotdVdPH(pH_dom)



    plt.show()
