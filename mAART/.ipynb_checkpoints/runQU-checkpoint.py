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

#start_time = time.perf_counter()
tL,QL,UL = calcQU.QULoops(spin,inclination,B_r,B_theta,B_phi,subkep_xi,beta_r,beta_phi,g_factor,r_start,phi_start)
#end_time = time.perf_counter()

# elapsed_time = end_time - start_time
# print(f"All necessary files calculated in: {elapsed_time:.2f} seconds\n")

print("Drawing QULoops")

plt.scatter(QL,UL,s=2)
plt.savefig("QU.png",dpi=200)
plt.close()

plt.scatter(tL,UL,s=2)
plt.savefig("UvT.png",dpi=200)
plt.close()

plt.scatter(tL,QL,s=2)
plt.savefig("QvT.png",dpi=200)
plt.close()

print("Success!")

save_dir = './QU_data/'
filename = f"i{inclination}_s{spin}_B({B_r},{B_theta},{B_phi})_subkep_{subkep_xi}_betar_{beta_r}_betaphi_{beta_phi}_rinit_{r_start}_phiinit_{phi_start}_gfactor_{g_factor}_QUdata.txt"
full_path = os.path.join(save_dir, filename)

np.savetxt(full_path, np.array([tL, QL, UL]).T)

'''
I0=calcQU.ImageQU(spin,inclination,B_r,B_theta,B_phi,subkep_xi,beta_r,beta_phi,g_factor,r_start,phi_start)

filename=f"i{inclination}_s{spin}_B({B_r},{B_theta},{B_phi})_betar_{beta_r}_betaphi_{beta_phi}_Image.h5"
h5f = h5py.File(os.path.join(save_dir, filename), 'w')

h5f.create_dataset('bghts0', data=I0)

h5f.close()

print("File ",filename," created.")
'''
