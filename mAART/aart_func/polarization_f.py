from aart_func import *
from params import *
from aart_func.intensity_f import gGas,gDisk#, lhat, Omegahat, Omegabar,Delta,Ehat,nuhat,urbar,uttilde


'''
HERE THE POLARIZATION IS COMPUTED FOR A GENERAL FLOW

'''
@jit()
def get_l_kepler(rad, a, sub_kep):
    return sub_kep * (rad**2 + a**2 - 2*a*np.sqrt(rad)) / (np.sqrt(rad)*(rad - 2) + a)

@jit()
def get_metrics(r, a):
    delta = r**2 - 2*r + a**2
    A_pot = (r**2 + a**2)**2 - a**2 * delta
    return delta, A_pot

@jit()
def get_p_photon(r,a,lamb,eta,redshift_sign,sqR,delta):
    pt = (r*(a**4 - 2*a*lamb*r + r**4 + a**2*(-delta + 2*r**2)))/(a**4*(-2 + r) + (-2 + r)*r**4 + a**2*(2*delta + r*(4 - delta + 2*(-2 + r)*r)))
    pr = (redshift_sign*sqR)/r**2
    pth = -(np.sqrt(eta)/r**2)
    pphi = ((2*a + lamb*(-2 + r))*r**2)/(a**4*(-2 + r) + (-2 + r)*r**4 + a**2*(2*delta + r*(4 - delta + 2*(-2 + r)*r)))
    return pt, pr, pth, pphi

@jit()
def get_boost(r,a,A_pot,delta,ut,ur,uphi):
    vr = (np.sqrt(A_pot)*ur)/(delta*ut)
    vphi = (np.sqrt(A_pot)*((np.sqrt(A_pot)*uphi)/r - (2*a*ut)/np.sqrt(A_pot)))/(np.sqrt(delta)*r*ut)
    return vr, vphi

@jit()
def get_pol(r,a,lamb,eta,redshift_sign,sqR,delta,A_pot,vr,vphi,ut,ur,uphi,Br,Bth,Bphi):
    B = np.sqrt(Br*Br + Bth*Bth + Bphi*Bphi)
    ft = (np.sqrt(A_pot)*vr*(-((Bphi*np.sqrt(eta))/r) - (2*a*Bth*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) + (np.sqrt(A_pot)*Bth*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (Bth*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) - (Bth*lamb*r*vr**2)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) - (Bth*lamb*r*vphi**2)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) - (Bth*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2))))/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2)) + (np.sqrt(A_pot)*vphi*((Br*np.sqrt(eta))/r + (2*a*Bth*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*Bth*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (Bth*redshift_sign*sqR*vphi**2)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) - (Bth*lamb*r*vphi*vr)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (Bth*lamb*r*vphi*vr)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) + (Bth*redshift_sign*sqR*vr**2)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2))))/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2))
    fr = (np.sqrt(delta)*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2))*(-((Bphi*np.sqrt(eta))/r) - (2*a*Bth*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) + (np.sqrt(A_pot)*Bth*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (Bth*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) - (Bth*lamb*r*vr**2)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) - (Bth*lamb*r*vphi**2)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) - (Bth*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2))))/(r*(vphi**2 + vr**2)*np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2)) + (np.sqrt(delta)*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2))*((Br*np.sqrt(eta))/r + (2*a*Bth*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*Bth*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (Bth*redshift_sign*sqR*vphi**2)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) - (Bth*lamb*r*vphi*vr)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (Bth*lamb*r*vphi*vr)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) + (Bth*redshift_sign*sqR*vr**2)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2))))/(r*(vphi**2 + vr**2)*np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2))
    fth = ((2*a*Br*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*Br*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) - (2*a*Bphi*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) + (np.sqrt(A_pot)*Bphi*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) - (Bphi*redshift_sign*sqR*vphi**2)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (Bphi*lamb*r*vphi*vr)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) - (Br*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (Br*lamb*r*vr**2)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (Br*lamb*r*vphi**2)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) - (Bphi*lamb*r*vphi*vr)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) + (Br*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) - (Bphi*redshift_sign*sqR*vr**2)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)))/(r*np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2))
    fphi = ((-((Bphi*np.sqrt(eta))/r) - (2*a*Bth*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) + (np.sqrt(A_pot)*Bth*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (Bth*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) - (Bth*lamb*r*vr**2)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) - (Bth*lamb*r*vphi**2)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) - (Bth*redshift_sign*sqR*vphi*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)))*((2*a*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) + (r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2))))/np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2) + (((Br*np.sqrt(eta))/r + (2*a*Bth*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*Bth*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (Bth*redshift_sign*sqR*vphi**2)/(np.sqrt(delta)*r*(vphi**2 + vr**2)) - (Bth*lamb*r*vphi*vr)/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (Bth*lamb*r*vphi*vr)/(np.sqrt(A_pot)*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)) + (Bth*redshift_sign*sqR*vr**2)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)*(vphi**2 + vr**2)))*((2*a*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) + (r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2))))/np.sqrt(eta/r**2 + ((2*a*lamb*vphi)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vphi)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (redshift_sign*sqR*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)) + (lamb*r*(vr**2 + vphi**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)))**2 + ((2*a*lamb*vr)/(np.sqrt(A_pot)*np.sqrt(delta)*np.sqrt(1 - vphi**2 - vr**2)) - (np.sqrt(A_pot)*vr)/(np.sqrt(delta)*r*np.sqrt(1 - vphi**2 - vr**2)) + (lamb*r*vphi*vr*(-1 + 1/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(A_pot)*(vphi**2 + vr**2)) + (redshift_sign*sqR*(vphi**2 + vr**2/np.sqrt(1 - vphi**2 - vr**2)))/(np.sqrt(delta)*r*(vphi**2 + vr**2)))**2)
    return ft/B, fr/B, fth/B, fphi/B


def KGas(r,thetad,a,lamb,eta,redshift_sign,sqR,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor):
    """
    Calculates the Penrose–Walker constant (Eq. 6 P1)
    """
    # 1. Metric at current r
    delta, A_pot = get_metrics(r, a)

    isco=rms(a)

    # 2. ISCO Constants (Conserved during plunge)
    L_isco = get_l_kepler(isco, a, sub_kep)
    delta_isco, A_pot_isco = get_metrics(isco, a)

    # Energy at ISCO
    denom_E = A_pot_isco/(isco**2) - (4*a*L_isco)/isco - (1 - 2/isco)*L_isco**2
    E_isco = np.sqrt(delta_isco / denom_E)

    # 3. Kinematics (Plunge Physics)
    omega_kepler = (a + (1 - 2/r)*(L_isco - a)) / (A_pot/(r**2) - (2*a*L_isco)/r)

    potential_term = A_pot/(r**2) - (4*a*L_isco)/r - (1 - 2/r)*L_isco**2 - delta/(E_isco**2)
    v_radial_phys = (r / delta) * np.sqrt(np.abs(potential_term))
    ur_hat = -(delta / r**2) * v_radial_phys * E_isco

    # Free-fall frame for mixing
    ur_bar = -np.sqrt(2*r*(r**2 + a**2)) / (r**2)
    omega_bar = (2*a*r) / A_pot

    # Mixing
    ur = ur_hat + (1 - betar) * (ur_bar - ur_hat)
    OT = omega_kepler + (1 - betaphi) * (omega_bar - omega_kepler)

    # 4. Normalization (uttilde)
    norm_num = 1 + (ur**2 * r**2) / delta
    norm_den = 1 - (r**2 + a**2) * OT**2 - (2/r) * (1 - a*OT)**2

    ut = np.sqrt(np.abs(norm_num / norm_den))
    uphi = ut * OT

    # r2 = r*r
    # r3 = r2*r
    # r4 = r2*r2
    # a2 = a*a
    # a3 = a2*a
    # a4 = a2*a2
    # delta2 = delta*delta
    # ut2 = ut*ut
    # ur2 = ur*ur
    # uphi2 = uphi*uphi

    pt, pr, pth, pphi = get_p_photon(r,a,lamb,eta,redshift_sign,sqR,delta)

    vr, vphi = get_boost(r,a,A_pot,delta,ut,ur,uphi)

    ft, fr, fth, fphi = get_pol(r,a,lamb,eta,redshift_sign,sqR,delta,A_pot,vr,vphi,ut,ur,uphi,Br,Bth,Bphi)

    k1 = r*(pt*fr - pr*ft + a*(pr*fphi - pphi*fr))

    k2 = -r*((r*r + a*a)*(pphi*fth - pth*fphi) - a*(pt*fth - pth*ft))

    return np.array((k1+1j*k2),dtype=np.complex128)


def KDisk(r,thetad,a,lamb,eta,redshift_sign,sqR,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor):
    """
    Calculates the Penrose–Walker constant (Eq. 6 P1)
    """
    # 1. Metric and Local Keplerian L
    delta, A_pot = get_metrics(r, a)
    L_local = get_l_kepler(r, a, sub_kep)

    # 2. Kinematics
    OT_kepler = (a + (1 - 2/r)*(L_local - a)) / (A_pot/(r**2) - (2*a*L_local)/r)
    ur_bar = -np.sqrt(2*r*(r**2 + a**2)) / (r**2)
    omega_bar = (2*a*r) / A_pot

    # Mixing
    ur = (1 - betar) * ur_bar
    OT = OT_kepler + (1 - betaphi) * (omega_bar - OT_kepler)

    # 3. Normalization
    norm_num = 1 + (ur**2 * r**2) / delta
    norm_den = 1 - (r**2 + a**2) * OT**2 - (2/r) * (1 - a*OT)**2

    ut = np.sqrt(np.abs(norm_num / norm_den))
    uphi = ut * OT

    pt, pr, pth, pphi = get_p_photon(r,a,lamb,eta,redshift_sign,sqR,delta)

    vr, vphi = get_boost(r,a,A_pot,delta,ut,ur,uphi)

    ft, fr, fth, fphi = get_pol(r,a,lamb,eta,redshift_sign,sqR,delta,A_pot,vr,vphi,ut,ur,uphi,Br,Bth,Bphi)

    k1 = r*(pt*fr - pr*ft + a*(pr*fphi - pphi*fr))

    k2 = -r*((r*r + a*a)*(pphi*fth - pth*fphi) - a*(pt*fth - pth*ft))

    return np.array((k1+1j*k2),dtype=np.complex128)


def kappa(a,thetao,Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor,n_band=0,save_to_disk=False):

    fnbands=path+"LensingBands_a_%f_i_%f.h5"%(a,thetao)

    h5f = h5py.File(fnbands,'r')
    grid = h5f[f'grid{n_band}'][:]
    mask = h5f[f'mask{n_band}'][:]
    h5f.close()

    fnrt=path+"Rays_a_%f_i_%f.h5"%(a,thetao)

    h5f = h5py.File(fnrt,'r')
    rs = h5f[f'rs{n_band}'][:]
    redshift_sign = h5f[f'sign{n_band}'][:]
    h5f.close()

    """
    Computes the linear polarization
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]

    #Conserved quantities
    lamb = -alpha*sin(thetao)
    eta = (alpha**2-a**2)*cos(thetao)**2+beta**2

    redshift_sign = redshift_sign[mask]

    Deltaaux=rs**2-2*rs+a**2

    sqR=np.sqrt((rs**2+a**2-a*lamb)**2-Deltaaux*(eta+(lamb-a)**2))

    r_p = 1+np.sqrt(1-a**2)
    isco=rms(a)

    #Redshift factor
    gTotal = np.zeros(rs.shape[0],dtype=np.complex128)

    gTotal[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco],Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor)
    gTotal[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco],Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor)

    gTotal[rs<=r_p] = 0

    gTotals = np.zeros(mask.shape,dtype=np.complex128)
    gTotals[mask] = gTotal

    polk = np.zeros(rs.shape[0],dtype=np.complex128)

    polk[rs>=isco]= KDisk(rs[rs>=isco],thetad,a,lamb[rs>=isco],eta[rs>=isco],redshift_sign[rs>=isco],sqR[rs>=isco],Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor)
    polk[rs<isco] = KGas(rs[rs<isco],thetad,a,lamb[rs<isco],eta[rs<isco],redshift_sign[rs<isco],sqR[rs<isco],Br,Bth,Bphi,betar,betaphi,sub_kep,gfactor)

    polk[rs<=r_p] = 0

    #polk = pw(rs, thetad, a, lamb, eta, redshift_sign, sqR, r_p)

    k1=np.real(polk)
    k2=np.imag(polk)

    #Electric vector polarization angle EVPA (Eq. 5 P1)
    nu=-(alpha+a*sin(thetao))

    #EVPA_d=sqrt((k1**2+k2**2)*(beta**2+nu**2))
    EVPA_d=(nu**2+beta**2)

    EVPA_i=(beta*k2-nu*k1)
    EVPA_j=(beta*k1+nu*k2)

    EVPA_i[rs<=r_p] = np.nan
    EVPA_j[rs<=r_p] = np.nan

    mask_d=EVPA_d>0
    EVPA_i[mask_d]/=EVPA_d[mask_d]
    EVPA_j[mask_d]/=EVPA_d[mask_d]

    PK = np.zeros(mask.shape,dtype=np.complex128) # ***Change: np.complex_ to np.complex128***
    EVPA_x = np.zeros(mask.shape)
    EVPA_y = np.zeros(mask.shape)

    EVPA_x[mask] = EVPA_i
    EVPA_y[mask] = EVPA_j
    PK[mask]     = polk

    '''
    plt.imshow((EVPA_x.real).reshape(6000,6000).T)
    plt.show()
    print(np.sum(EVPA_x.real))
    '''

    data_out = {
        # "isco": isco,
        # "PK": PK,
        "EVPA_x": EVPA_x,
        "EVPA_y": EVPA_y,
        "gTotal": gTotals
    }

    if save_to_disk:
        filename=path+"Polarization_a_%f_i_%f_Br_%f_Bth_%f_Bphi_%f_betar_%f_betaphi_%f_sub_kep_%f.h5"%(np.round(a,6),np.round(thetao,6),np.round(Br,6),np.round(Bth,6),np.round(Bphi,6),np.round(betar,6),np.round(betaphi,6),np.round(sub_kep,6))

        # h5f = h5py.File(filename, 'a')
        with h5py.File(filename, 'a') as h5f:
            if f'gTotal_n{n_band}' not in h5f:
                #h5f.create_dataset('isco', data=isco)
                h5f.create_dataset(f'gTotal_n{n_band}', data=gTotals)
                #h5f.create_dataset('PK',     data=PK)
                h5f.create_dataset(f'EVPA_x_n{n_band}', data=EVPA_x)
                h5f.create_dataset(f'EVPA_y_n{n_band}', data=EVPA_y)


                if n_band==0:
                    print("File ",filename," created.")
                else:
                    print("File ",filename,f" is being updated with band n = {n_band} values.")
            else:
                print(f"n = {n_band} data already in", filename,".")

        h5f.close()
        return 0

    # OTHERWISE, JUST RETURN THE DATA
    return data_out

