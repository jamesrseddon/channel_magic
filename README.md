
# channel_magic v2.0


    Channel Magic provides a set of tools for calculating measures of magic
    for multi-qubit channels, for use with the CVX convex optimisation 
    package in MATLAB.
    Copyright (C) 2019-2022  James Seddon 
    
    MIT License

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
    


Email: jseddonquantum@gmail.com



This code accompanies the article
[Seddon2019]
"Quantifying magic for multi-qubit operations",
by James R. Seddon and Earl T. Campbell,
Proc. R. Soc. A. 475:20190251 (2019)
https://doi.org/10.1098/rspa.2019.0251
https://arxiv.org/abs/1901.03322

as well as the thesis
[Seddon2021]
"Advancing classical simulators by measuring the magic of quantum computation", 
by James R. Seddon,
University College London (2021).
This thesis will shortly be available via UCL Open Access, and this readme
will be updated when a link is available.

The code also allows computation of several of the measures of magic 
defined in
[Seddon2021B]
"Quantify quantum speedups: improved classical simulation from tighter magic monotones",
by James R. Seddon, Bartosz Regula, Hakop Pashayan, Yingkai Ouyang and Earl T. Campbell,
PRX Quantum 2, 010345 (2021)
https://doi.org/10.1103/PRXQuantum.2.010345
https://arxiv.org/abs/2002.06181

The MATLAB tools released in v1.0 of this repository were primarily aimed
at calculating the channel robustness and magic-increasing capacity, as 
described in Seddon and Campbell (2019), for two-qubit CPTP maps in the most 
general case, and for up to five qubits in the case of diagonal channels.
The repository also contains several data files Copyright (C) 2017 Mark Howard, see CITATION.md for details.


v2.0 includes code to compute or bound further measures of magic for
channels defined in the thesis cited above [Seddon2021], as well as measures
of magic for states discussed in [Seddon2021B].


------------------------------
# Prerequisites and Installation

Code is written in MATLAB, and requires the free convex optimisation
package CVX.
Downloads and documentation for CVX can be found at cvxr.com.

The repository also contains two Mathematica notebooks. These relate to
the plotting of specific figures in the thesis [Seddon2021], and
Mathematica is not needed for using main functionality.

To install, simply clone this repository and add to your MATLAB path.

This code has been tested using:
MATLAB version 9.0 (R2016a) on Windows, with CVX version 2.1, Build 1116.
MATLAB version 9.5 (R2018b) on Windows, with CVX version 2.1, Build 1127.


For full functionality the following freely available third-party code must
be in the user's MATLAB path:
(1) From Bartosz Regula's code accompanying PRX Quantum 2, 010345 (2021)
    https://bartoszregula.me/code/magic, the functions:
    'Lambda' for computing dyadic negativity, and
    'RobMagG' for computing generalised robustness. 
(2) The following MATLAB functions by Toby Cubitt, available at
    https://www.dr-qubit.org/matlab.html
    'randRho','randU' and 'syspermute'.


------------------------------
# Working examples

The root directory contains several stand-alone scripts that illustrate 
the basic functionality of the code:

channel_robustness.m
--------------------
Calculates the channel robustness for several simple single-qubit maps.

magic_capacity.m
----------------
Calculates magic capacity for several examples.

diagonal_channels.m
-------------------
Calculates both magic measures for some example diagonal operations on up 
to 5 qubits.

load_robustness_files.m
-----------------------
A number of the functions in this repository involve linear systems of the 
form Ax = b. This script loads the data files containing the necessary 
A matrices into the workspace. It also loads cell arrays containing affine 
space representations of stabiliser states, which are needed for 
calculating capacity in the case of diagonal channels.

-------------------------------
# Reproducing numerical results

[Seddon2019]

The 'scripts' directory contains scripts for reproducing the numerical 
results presented in arXiv:1901.03322.

The results for single-qubit rotations subject to amplitude damping can be 
reproduced using '/scripts/noisy_channel_sweep.m'.

Results for multi-control phase gates and random diagonal gates can be 
reproduced using '/scripts/multicontrol_phase_gates.m' and 
'/scripts/random_phase_searches.m' respectively.

Results for tensor product Z-rotations can be calculated using 
'/scripts/Z_rotation_sweep.m'.

[Seddon2021]

The subdirectory 'scripts/thesis_6_2' includes MATLAB scripts for 
generating the data used to plot figures 6.8 to 6.10 from Section 6.2, 
Chapter 6 of the PhD thesis cited above [Seddon2021].

the directory 'mathematica_notebooks' contains notebooks that can be used
to reproduced figures 6.6, 6.7 and 6.10 of the same chapter.






