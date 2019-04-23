% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
%  Copyright © 2019 James R. Seddon
%  This file is part of project: Channel Magic
%
%%% Load affine space data and A matrices for linear optimisation %%%

%%%% Load A matrices for linear system %%%%
load('Amat1.mat'); % Load A matrix for all 1-qubit stabiliser states.
A_mat_1 = Amat1;
clear A;
load('Amat2.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_2 = A;
clear A;
load('Amat3.mat'); % Load A matrix for all 3-qubit stabiliser states.
A_mat_3 = A;
clear A;
load('Amat4.mat'); % Load A matrix for all 4-qubit stabiliser states.
A_mat_4 = A;
clear A;
load('Amat5.mat'); % Load A matrix for all 5-qubit stabiliser states.
A_mat_5 = A; % N.B. This matrix takes up ~1.2 GB memory when loaded.
clear A;

load('Amat1_Clifford_Pauli_reset.mat') % Load A matrix for all Choi
% matrices corresponding to single-qubit Cliffords and Pauli reset
% channels (CPR).

%%% Load affine space data %%%
load('affine_spaces_one_qubit.mat');
load('affine_spaces_two_qubit.mat');
load('affine_spaces_three_qubit.mat');
load('affine_spaces_four_qubit.mat');
load('affine_spaces_five_qubit.mat');
load('affine_space_HEADERS.mat'); % See this array for header text for 
                                  % for affine space data arrays.
                                  
  