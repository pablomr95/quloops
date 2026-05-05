from aart_func import *

print("\nThanks for using mAART")
print("Copyright (C) 2026, P. Ruales, A. Cardenas-Avendano, D. Gates\n")

#--------------------------------------#
# Magnetic field in spherical polar coordinates in the local frame

Br=0
Bth=1
Bphi=0

#--------------------------------------#

#BH's Spin
spin_case=0.94
#Observer's inclination
i_case=20

#Velocity Profile for the gas

#Sub-Kepleniarity param
sub_kep=0.95
#Radial velocity param
betar=0.98
#Angular velocity param
betaphi=0.98

#Initial conditions for the trajectory
r_initial = 11
phi_initial = 0

#The power of the redshift factor
gfactor=2

#isco = rms(spin_case)

'''
MIT license
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.
'''
