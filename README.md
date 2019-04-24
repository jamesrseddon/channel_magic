
# channel_magic v1.0


    Channel Magic provides a set of tools for calculating measures of magic for multi-qubit channels, 
    for use with the CVX convex optimisation package in MATLAB.
    Copyright (C) 2019  James Seddon 
    
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
    
    


Email: jamesseddon@physics.org



This code accompanies the article "Quantifying magic for multi-qubit operations", by James R. Seddon and Earl T. Campbell, 2019.

A pre-print version of the article can be found here: https://arxiv.org/abs/1901.03322

These MATLAB tools are primarily aimed at calculating the channel robustness and magic-increasing capacity, as described in the above paper, for two-qubit CPTP maps in the most general case, and for up to five-qubit diagonal channels.

The repository also contains several data files Copyright (C) 2017 Mark Howard, see CITATION.md for details.

------------------------------
# Prerequisites and Installation

Code is written in MATLAB, and requires the free convex optimisation package CVX.
Downloads and documentation for CVX can be found at cvxr.com.

This code has been tested using:
MATLAB version 9.0 (R2016a) on Windows, with CVX version 2.1, Build 1116.
MATLAB version 9.5 (R2018b) on Windows, with CVX version 2.1, Build 1127.

To install, simply clone this repository and add to your MATLAB path.

------------------------------
# Working examples

The root directory contains several stand-alone scripts that illustrate the basic functionality of the code:

channel_robustness.m
--------------------
Calculates the channel robustness for several simple single-qubit maps.

magic_capacity.m
----------------
Calculates magic capacity for several examples.

diagonal_channels.m
-------------------
Calculates both magic measures for some example diagonal operations on up to 5 qubits.

load_robustness_files.m
-----------------------
A number of the functions in this repository involve linear systems of the form Ax = b. This script loads the data files containing the necessary A matrices into the workspace. It also loads cell arrays containing affine space representations of stabiliser states, which are needed for calculating capacity in the case of diagonal channels.

-------------------------------
# Reproducing numerical results

The 'scripts' directory contains scripts for reproducing the numerical results presented in arXiv:1901.03322.

The results for single-qubit rotations subject to amplitude damping can be reproduced using '/scripts/noisy_channel_sweep.m'.

Results for multi-control phase gates and random diagonal gates can be reproduced using '/scripts/multicontrol_phase_gates.m' and '/scripts/random_phase_searches.m' respectively.

Results for tensor product Z-rotations can be calculated using '/scripts/Z_rotation_sweep.m'.

