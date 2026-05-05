from aart_func import *
from params import *
from paramsQU import *

# This function runs lensingbands.py, raytracing.py, polarization.py and calculates the Q and U Stokes linear polarization for each photon

def QULoops(spin_case,thetao,Br,Bth,Bphi,sub_kep,betar,betaphi,gfactor,r_initial,phi_initial,n_turns=3,t_max=2000):

    traj = polp.generate_trajectory(spin_case,betar,betaphi,sub_kep,r_initial,phi_initial,n_turns,t_max)

    pol_data = polarizationk.kappa(spin_case,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor)

    evpax = pol_data['EVPA_x']
    evpay = pol_data['EVPA_y']
    gTotal = pol_data['gTotal']

    t_n, Q_n0, U_n0, t_ray_n0 = polp.process_polarization(thetao, spin_case, Br, Bth, Bphi, betar, betaphi, sub_kep, gfactor, traj, evpax, evpay, gTotal)

    return t_n, traj[1], traj[2], Q_n0, U_n0, t_ray_n0
