from aart_func import*
from params import*
from paramsQU import*
import time

spin = np.round(spin_case,6)
inclination = np.round(i_case,6)
B_r = np.round(Br,6)
B_theta = np.round(Bth,6)
B_phi = np.round(Bphi,6)
subkep_xi = np.round(sub_kep,6)
beta_r = np.round(betar,6)
beta_phi = np.round(betaphi,6)
g_factor = np.round(gfactor,6)
r_start = np.round(r_initial,6)
phi_start = np.round(phi_initial,6)


I0=calcQU.ImageQU(spin,inclination,B_r,B_theta,B_phi,subkep_xi,beta_r,beta_phi,g_factor,r_start,phi_start,n_band=0)
I1=calcQU.ImageQU(spin,inclination,B_r,B_theta,B_phi,subkep_xi,beta_r,beta_phi,g_factor,r_start,phi_start,n_band=1)
I2=calcQU.ImageQU(spin,inclination,B_r,B_theta,B_phi,subkep_xi,beta_r,beta_phi,g_factor,r_start,phi_start,n_band=2)

save_dir = './QU_data/'
filename=f"i{inclination}_s{spin}_B({B_r},{B_theta},{B_phi})_betar_{beta_r}_betaphi_{beta_phi}_Image.h5"
h5f = h5py.File(os.path.join(save_dir, filename), 'w')

h5f.create_dataset('bghts0', data=I0)
h5f.create_dataset('bghts1', data=I1)
h5f.create_dataset('bghts2', data=I2)

h5f.close()

print("File ",filename," created.")
