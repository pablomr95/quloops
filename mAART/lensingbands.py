from aart_func import *
from paramsQU import *

thetao=np.round(i_case*np.pi/180,6)

print("Computing the lensing bands")
lb.lb(spin_case,thetao)
