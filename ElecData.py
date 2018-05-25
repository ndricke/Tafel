import numpy as np


"""
ElecData is a class that generates electrochemical data for fitting with proposed rate laws 
"""
class ElecData(object):
    
    def __init__(s, pH_points, pH_range, V_points, V_range):
        s.pH_dom = np.linspace(pH_range[0],pH_range[1],pH_points)
        s.V_dom = np.linspace(V_range[0],V_range[1],V_points)
    
    @staticmethod
    def linearPH( V_slope, V_int, pH_dom):
        onsetV = -1.*V_slope*pH_dom - V_int
        return onsetV

    @staticmethod
    def linearTafel(logI_slope, logI_int, V_dom):
        tafel_curve = -1.*logI_slope*V_dom - logI_int
        return tafel_curve

