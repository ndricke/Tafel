import numpy as np
import matplotlib.pyplot as plt

#prev = partially reversible cycle, where every step is reversible except for one
#this function assumes that step 4 is irreversible
def prevCycle4(k1,k2,k3,k4,kn1,kn2,kn3):
    num = k1*k2*k3
    denom = kn1*kn2*kn3+ kn1*kn2*k4+ kn1*k3*k4+ k2*k3*k4+ kn2*kn3*k1+ kn2*k1*k4+ k1*k3*k4+ \
            kn3*k1*k2+ k1*k2*k4+ k1*k2*k3
    c4 = denom/num
    r4 = k4*num/denom
    c3 = c4*(kn3+k4)/k3
    c2 = (kn2*c3+k4*c4)/k2
    c1 = (kn1*c2+k4*c4)/k1
    return r4,c1,c2,c3,c4

#should be given 2n-1 rate constants, where n is the number of steps in the cycle
def genPrevCycle(k)
    pass

def prev4pOH(OH,ks,OHdep):
    for i in range(len(ks)):
        ks[i] *= OH**OHdep[i]
    return prevCycle4(*ks)

#Want to scan over pOH dependent steps for log(j) dependence on current
#Let's assume that all of the reversible steps are pOH dependent in the reverse direction

OHdep4 = [0,0,0,0,1,1,1]
ks = [1.,1.,1.,1.,1.,1.,1.]

OH = 0.2




