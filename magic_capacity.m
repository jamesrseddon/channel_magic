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
%%%% Script to calculate magic capacity for example channels %%%%

%%% Define some operations %%%

Tgate = [ 1 0;
               0 exp(1i*pi/4)];
           
Hadamard = (1/sqrt(2))*[1 1;
                        1 -1];
                    
% Amplitude damping noise channel with noise parameter p, Kraus operators:
%   K_1 = |0><0| + sqrt(1-p)|1><1|
%   K_2 = sqrt(p)|0><1|

amp_damp_0pt1 = amplitudeDamping(0.1);

% Amplitude damping with noise p followed by X-rotation by angle theta. 
params.p = 0.1;
params.theta = pi/8;
amp_then_Xrotate = constructChannel('AMP_X',params);


%%%% Load A matrices for linear system %%%%
load('Amat2.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_2 = A;
clear A;

[capacity_amp, result_table_amp] = findCapacity(amp_damp_0pt1,...
                                        A_mat_2,A_mat_2,'low','',...
                                        'capacity_results_amp_damp.csv');
                                    
[capacity_ampX, result_table_ampX] = findCapacity(amp_then_Xrotate,...
                                        A_mat_2,A_mat_2,'low','',...
                                        'capacity_results_ampX.csv');
                                    
result_array = {'Operation','Magic capacity';
                'Amplitude damping channel',capacity_amp;
                'Amplitude damped X-rotation',capacity_ampX}

                                    