from aart_func import *
#from params import *
from paramsQU import *

thetao=np.round(i_case*np.pi/180,6)

polarizationk.kappa(spin_case,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor,n_band=0,save_to_disk=True)
polarizationk.kappa(spin_case,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor,n_band=1,save_to_disk=True)
polarizationk.kappa(spin_case,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor,n_band=2,save_to_disk=True)
