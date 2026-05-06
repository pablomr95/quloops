[![2602.09102](https://img.shields.io/badge/arXiv-2602.09102-b31b1b.svg)](https://arxiv.org/abs/2602.09102) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/pablomr95/quloops/blob/main/License.txt)

# pAART #

"pAART" is a version of AART (https://github.com/iAART/aart), modified to calculate the Stokes Q and U parameters along the trajectory that an orbiting emitter takes within the disk in the equatorial plane. It is ready to calculate and visualize the linear polarization values in the observer's screen. Below is an excerpt from AART's README file.

# Adaptive Analytical Ray Tracing (AART) #

AART is a numerical framework that exploits the integrability properties of the Kerr spacetime to compute high-resolution black hole images and their visibility amplitude on long interferometric baselines. It implements a non-uniform adaptive grid on the image plane suitable to study black hole photon rings (narrow ring-shaped features predicted by general relativity but not yet observed). 

The code, described in detail in Ref. [2], implements all the relevant equations required to compute the appearance of equatorial sources on the (far) observer's screen. We refer the Reader to Refs. [3-5] for the derivations and further details. Through the code, the equations are mentioned as Pi Eq. N, which means Eq. N in Ref. [i]. 

The use of AART in scientific publications must be properly acknowledged. Please cite:

_______
Cardenas-Avendano, A., Lupsasca, A. & Zhu, H. "Adaptive Analytical Ray Tracing of Black Hole Photon Rings." Phys. Rev. D 107, 043030, 2023. [arXiv:2211.07469](https://arxiv.org/abs/2211.07469) 

and

P. Ruales, D. E. A. Gates, and A. Cárdenas-Avendaño. Polarization signatures of inspiraling hotspots around Kerr Black Holes. Phys. Rev. D (Accepted May, 2026). [arXiv:2602.09102] (https://arxiv.org/abs/2602.09102)
_______

We also request that AART modifications or extensions leading to a scientific publication be made public as free software. 

<center> <em>Feel free to use images and movies produced with this code (with attribution) for your next presentation! </em> </center>

_______
![GitHub last commit](https://img.shields.io/github/last-commit/pablomr95/quloops)
_______

## AART's Components ##

* **Lensing Bands**: The main functions are located in <em>lb_f.py</em> : This module computes the Bardeen's coordinates inside the so-called lensing bands (currently it only computes ($0\le n\le 2$), and the extension to a higher n is possible: just compy the structure of the code and add the desired n number) on a Cartesian grid with different resolutions. 

* **Analytical Ray-Tracing**: The main functions are located in  <em>raytracing_f</em>: For a given location in the Bardeen's plane ($\alpha,\beta$), it computes where it lands in the equatorial plane ($t,r,\theta=\pi/2,\phi$) in Boyer-Lindquist coordinates. The implementatio does it per lensing band. 

## Dependencies ##

#### Python Libraries: 

All the dependencies are located in the <em>init.py</em> file. Most of the libraries will come natively with anaconda (e.g., numpy, scipy >=1.8, matplotlib, multiprocessing, skimage) but some may not. 

To install all requirements*, run

 <code> pip install -r requirements.txt </code>

or, if using anaconda,

 <code> conda install --yes --file requirements.txt </code>

You can also install any missing packages by running
  
<code> pip install "package_name" </code>
  
or, if using anaconda, search for the missing packages and run, e.g. for h5py (Read and write HDF5 files from Python,) 
  
<code> conda install -c anaconda h5py</code>

Sometimes scipy does not update automatically to the latest version. If that is the case, you may want to type 

<code> pip install -U scipy</code>

Some users have experienced an issue with <em>imageio.v2</em>, as it is not found. To solve this issue please type:

<code> python -m pip install --upgrade pip </code>

<code> pip install imageio --upgrade </code>

<em>*Thanks to @prestonyun for suggesting this simplification.</em> 

## How to run AART and calculate QU loops ##

### As a python package:

Simply [pip](https://pypi.org/project/aart/) install it like this:

<code> pip install aart </code>

Then follow the notebook: 

<em>Example_Notebook_InspiralingHotspot.ipynb</em>

To change the paramaters modify the file and <em>paramsQU.py</em>.


#### Lensing Bands: 

The lensing bands are computed by simply running

  <code> python lensingbands.py </code>
  
The result will be stored in a HDF5 file that contains the values of the Bardeen's coordinates within each lensing band. The datasets inside the resulting file are:

* alpha: The coordinate alpha of the critical curve. The parameter <em>npointsS</em> controls the number of points used for the computation of the critical curve)
* beta: The coordinate beta of the critical curve. 

* hull\_ni: The points for the inner convex hull of the nth band. Note that hull\_0i corresponds to the location of the apparent horizon. 
hull\_ne: The points for the outer convex hull of the nth band. Note that hull_0e corresponds to edges of the domain.  
* gridn: The point within the nth lensing band. 
* Nn: Number of points within the nth lesing band.
* limn: The grids are cartesian and symmetric around zero. This data sets tells the limits of the grid. 


#### Ray Tracing: 

To compute the equitorial radius, angle, and emission time of a photon, we perform a backward ray-tracing from the observer plane. By running the following, we evaluate the source radius, angle, and time within the grid from each lensing bands:

  <code> python raytracing.py </code>
  
The result will be stored in a HDF5 file that contains source radius, angle, time, as well as the radial component of the four momentum at the equitorial plane, for lensing bands n=0,1,2. The datasets inside the resulting file are:

* rsn: The value of the r Boyer-Lindquist coordinate for the nth lensing band. It follows the order of the lensing band. 
* tn: The value of the t Boyer-Lindquist coordinate for the nth lensing band. It follows the order of the lensing band. 
* phin: The value of the \phi Boyer-Lindquist coordinate for the nth lensing band. It follows the order of the lensing band. 
* signn: The sign of the radial momentum of the emitted photon in the nth lensing band. 



#### Polarization:

The linear polarization of a given configuration of the magnetic field can be computed by

  <code> python polarization.py </code>
  
  The resulting datasets inside the resulting file are:

 
  * EVPA_x: The x-component of the the electric-vector position angle.
  * EVPA_y: The y-component of the the electric-vector position angle.

## Limitations and known possible performance bottlenecks ##

* This code has only been tested on Mac OS (M1 and Intel) and on Ubuntu. 

* If you want to run a retrograde disk, you will have to apply symmetry arguments. In other words, run the positive spin case ($-a$), flip the resulting lensing bands and rays, and then compute the intensity on each pixel. Note that the characteristic radii need also to be modified. We plan to add this feature in a future version. 

* The Radon cut does not smoothly goes to zero. This is sometimes clear from the visamp, where you can see an extra periodicity (wiggle) on each local maxima. To solve this issue, increase the FOV of the $n=0$ image by providing a larger value for the variable <em>limits</em> in <em>params.py</em>. You can also modify the percentage of points used in <em>npointsfit</em> in <em>visamp_f.py</em>.

* Producing the lensing bands is taking too long. Sometimes, in particular for larger inclination values, computing the contours of the lensing bands and the points within it, takes a long time. The calculation can be made faster, but less accurate if you decrease the number of points used to compute the contours, i.e., by decreasing the value of the variable <em>npointsS</em> in <em>params.py</em>. It is faster to compute the convex Hull instead of the concave Hull (alpha shape), but then you will have to check that your are not missing points (having extra points is not an issue with the analytical formulae, as the results are masked out). If using the convex is okay, then you can also change the function <em>in_hull</em> in <em>lb_f.py</em> to use <em>hull.find_simplex</em> instead of <em>contains_points</em>.


## Authors ##

### Current Developers ###

- Pablo Ruales (rualpm25 [at] wfu [dot] edu)
- Alejandro Cardenas-Avendano

### Former Developers ###
- Lennox Keeble
- Hengrui Zhu
- Alex Lupsasca


## References ##

[1] P. Ruales, D. E. A. Gates, and A. Cárdenas-Avendaño. Polarization signatures of inspiraling hotspots around Kerr Black Holes. Phys. Rev. D (Accepted May, 2026). [arXiv:2602.09102] (https://arxiv.org/abs/2602.09102)

[2] Cardenas-Avendano, A., Lupsasca, A. & Zhu, H. Adaptive Analytical Ray Tracing of Black Hole Photon Rings. Physical Review D, 107, 043030, 2023. [arXiv:2211.07469](https://arxiv.org/abs/2211.07469)

[3] Gralla, S. E., & Lupsasca, A. (2020). Lensing by Kerr black holes. Physical Review D, 101, 044031.

[4] Gralla, S. E., & Lupsasca, A. (2020). Null geodesics of the Kerr exterior. Physical Review D, 101, 044032.

[5] Gralla, S. E., Lupsasca, A., & Marrone, D. P. (2020). The shape of the black hole photon ring: A precise test of strong-field general relativity. Physical Review D, 102, 124004.

## MIT License

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
