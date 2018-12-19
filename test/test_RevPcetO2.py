import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel






dom = 4000
pH_dom = np.linspace(3,10,dom)
V_dom = np.linspace(0.,0.1,dom)

cyc1 = ElecMech.RevPcetO2((1.,1.,1.,1.), 0.)
cyc2 = ElecMech.RevPcetO2((1.,1.,1.,1.), 0.)

sim = SimTafel.SimTafel(pcet2)
