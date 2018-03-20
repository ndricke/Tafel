import numpy


class SemiCyc3(object):
    f = 38.684
    ep = 10.**-3

    def __init__(s,k1,k3,kn3,P_O2,pH,V):
        s.k1 = k1
        s.k3 = k3
        s.kn3 = kn3
        s.P_O2 = P_O2
        s.setV(V) #we want to make sure V and expV are always the same
        s.setpH(pH) #likewise for [H+] and pH

    def setV(s,V):
        s.V = V
        s.expV = np.exp(-V*f)

    def setpH(s,pH):
        s.pH = pH
        s.H = 10.**(-pH)

    def rate(s):
        num = s.k1*s.k3*s.P_O2*s.H*s.expV
        denom = k1*s.P_O2 + k3*s.H
        return num/denom



class MechKinetics(object):

    def __init__(s,mech, ep = 10**-3)
        s.mech = mech #An already initialized kinetic mechanism class
        s.ep = ep

    def varV(s,V):
        s.mech.setV(V)
        return s.mech.rate() - s.ep

    def varpH(s,pH):
        s.mech.setpH(pH)
        return s.mech.rate() - s.ep

    def varO2(s,P_O2):
        s.mech.P_O2 = P_O2
        return s.mech.rate() - s.ep

    #Scans over a concentration of pH, O2, ect (using varpH or varO2 to set)
    #Finds potential measured where current is measured at close enough to 0
    def concenScan(s,conc_range, conc_func):
    V_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
    for i in H:
        cyc.H = i
        V_sol = opt.root(cyc.varV,V_list[-1],method='hybr')
        V_list.append(V_sol.x[0])
    V_list.pop(0)
    gradV = np.gradient(V_list,h[1]-h[0])




if __name__ == "__main__":

    


    V = 0.
    pH_range = [-2.,8]
    h = np.linspace(pH_range[0],pH_range[1],200)
    H = 10.**(-1.*h)
    OH = 10.**(h-14.)
    V_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
    for i in H:
        cyc.H = i
        V_sol = opt.root(cyc.varV,V_list[-1],method='hybr')
        V_list.append(V_sol.x[0])
    V_list.pop(0)
    gradV = np.gradient(V_list,h[1]-h[0])
    print gradV
