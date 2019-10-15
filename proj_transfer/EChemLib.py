import pandas as pd
import numpy as np



def convertHE(V, pH, H_electrode):
    """
    Converts SHE to RHE and vice versa
    E(RHE) = E(SHE) + 0.059pH
    parameters:
    V (numpy 1d array) input voltage
    pH (numpy 1d array) pH at each voltage
    H_electrode (string) desired hydrogen electrode, either RHE or SHE

    returns:
    V_conv (numpy 1d array)
    """

    conv_sign = {"SHE":1.0, "RHE":-1.0}
    V_conv = V - 0.059*pH*conv_sign[H_electrode]
    return V_conv
