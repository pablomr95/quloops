from aart_func import*
from params import*
from paramsQU import*

spin = np.round(spin_case,6)
inclination = np.round(i_case*np.pi/180,6)
B_r = np.round(Br,6)
B_theta = np.round(Bth,6)
B_phi = np.round(Bphi,6)
subkep_xi = np.round(sub_kep,6)
beta_r = np.round(betar,6)
beta_phi = np.round(betaphi,6)
g_factor = np.round(gfactor,6)
r_start = np.round(r_initial,6)
phi_start = np.round(phi_initial,6)


t_traj, r_traj, phi_traj, Q, U, t_ray = calcQU.QULoops(spin,inclination,B_r,B_theta,B_phi,subkep_xi,beta_r,beta_phi,g_factor,r_start,phi_start,n_turns=N_turns,t_max=T_max)

save_dir = './QU_data/'
filename = f"QUdata_i{i_case}_s{spin}_B({B_r},{B_theta},{B_phi})_subkep_{subkep_xi}_betar_{beta_r}_betaphi_{beta_phi}_rinit_{r_start}_phiinit_{phi_start}_gfactor_{g_factor}.h5"
full_path = os.path.join(save_dir, filename)

h5f = h5py.File(full_path, 'w')

h5f.create_dataset('t_traj', data=t_traj)
h5f.create_dataset('r_traj', data=r_traj)
h5f.create_dataset('phi_traj', data=phi_traj)
h5f.create_dataset('Q', data=Q)
h5f.create_dataset('U', data=U)
h5f.create_dataset('t_ray', data=t_ray)

h5f.close()

print("File ",full_path," created.")

