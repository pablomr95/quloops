from aart_func import *
from params import *
# from paramsQU import *

# This function runs lensingbands.py, raytracing.py, polarization.py and calculates the Q and U Stokes linear polarization for each photon

def QULoops(spin_case,i_case,Br,Bth,Bphi,sub_kep,betar,betaphi,gfactor,r_initial,phi_initial):

    thetao=np.round(i_case*np.pi/180,6)

    # print("Computing the lensing bands")
    # lb.lb(spin_case,thetao)

    # print("Ray-tracing")

    fnbands=path+"LensingBands_a_%f_i_%f.h5"%(spin_case,thetao)

    print("Reading file: ",fnbands)

    h5f = h5py.File(fnbands,'r')

    #N0=int(h5f["N0"][0])
    supergrid0=h5f['grid0'][:]
    mask0=h5f['mask0'][:]
    supergrid1=h5f['grid1'][:]
    mask1=h5f['mask1'][:]
    supergrid2=h5f['grid2'][:]
    mask2=h5f['mask2'][:]
    h5f.close()

    # rt.rt(supergrid0,mask0,supergrid1,mask1,supergrid2,mask2,spin_case,thetao)
    # print("A total of",supergrid0.shape[0]+supergrid1.shape[0]+supergrid2.shape[0],"photons were ray-traced")

    fnrt=path+"Rays_a_%f_i_%f.h5"%(spin_case,thetao)

    print("Reading file: ",fnrt)

    h5f = h5py.File(fnrt,'r')

    rs0=h5f['rs0'][:]
    phi0 = h5f['phi0'][:]
    sign0=h5f['sign0'][:]
    h5f.close()

    pol_data = polarizationk.kappa(supergrid0,mask0,rs0,sign0,spin_case,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor,save_to_disk=False)

    evpax = pol_data['EVPA_x']
    evpay = pol_data['EVPA_y']
    gTotal = pol_data['gTotal']

    traj = polp.generate_trajectory(spin_case,betar,betaphi,sub_kep,r_initial,phi_initial)

    return polp.process_polarization(thetao, spin_case, Br, Bth, Bphi, betar, betaphi, sub_kep,gfactor,traj,rs0,phi0,evpax,evpay,gTotal)

def ImageQU(spin_case,i_case,Br,Bth,Bphi,sub_kep,betar,betaphi,gfactor,r_initial,phi_initial):

    thetao=np.round(i_case*np.pi/180,6)

    # print("Computing the lensing bands")
    # lb.lb(spin_case,thetao)

    # print("Ray-tracing")

    fnbands=path+"LensingBands_a_%f_i_%f.h5"%(spin_case,thetao)

    print("Reading file: ",fnbands)

    h5f = h5py.File(fnbands,'r')

    #N0=int(h5f["N0"][0])
    supergrid0=h5f['grid0'][:]
    mask0=h5f['mask0'][:]
    supergrid1=h5f['grid1'][:]
    mask1=h5f['mask1'][:]
    supergrid2=h5f['grid2'][:]
    mask2=h5f['mask2'][:]
    h5f.close()

    # rt.rt(supergrid0,mask0,supergrid1,mask1,supergrid2,mask2,spin_case,thetao)
    # print("A total of",supergrid0.shape[0]+supergrid1.shape[0]+supergrid2.shape[0],"photons were ray-traced")

    fnrt=path+"Rays_a_%f_i_%f.h5"%(spin_case,thetao)

    print("Reading file: ",fnrt)

    h5f = h5py.File(fnrt,'r')

    rs0=h5f['rs0'][:]
    sign0=h5f['sign0'][:]
    h5f.close()

    pol_data = polarizationk.kappa(supergrid0,mask0,rs0,sign0,spin_case,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor,save_to_disk=False)

    gTotal = pol_data['gTotal']

    #rH = 1 + np.sqrt(1 - spin_case**2)

    #gTotal[rs0==np.nan] = 0.0

    # SU's parameters for the envelope
    # Just used for the analytical profiles
    gammap=-3/2
    mup=1-sqrt(1-spin_case**2)
    sigmap=1/2

    return gTotal**3.0*ilp.profile(rs0,spin_case,gammap,mup,sigmap)
