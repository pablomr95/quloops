from aart_func import *
from params import *
from paramsQU import *

print("Ray-tracing")

thetao=np.round(i_case*np.pi/180,6)

fnbands=path+"LensingBands_a_%f_i_%f.h5"%(spin_case,thetao)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

supergrid0=h5f['grid0'][:]
mask0=h5f['mask0'][:]
supergrid1=h5f['grid1'][:]
mask1=h5f['mask1'][:]
supergrid2=h5f['grid2'][:]
mask2=h5f['mask2'][:]
h5f.close()

rt.rt(supergrid0,mask0,supergrid1,mask1,supergrid2,mask2,spin_case,thetao)
print("A total of",supergrid0.shape[0]+supergrid1.shape[0]+supergrid2.shape[0],"photons were ray-traced")

