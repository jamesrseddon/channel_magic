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
%%%% script for noisy channel sweep %%%%%

%%%% INPUTS %%%%

%%%% Load A matrices %%%
load('Amat2.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_2 = A;
clear A;

load('Amat1_Clifford_Pauli_reset.mat'); % Load A matrix for CPR Choi states.

A_matrices = {A_mat_2; A_mat_CP_1}; % List of A matrices. In configuration
                                    % set below, A matrices are indexed by
                                    % their position in this list.
                                    
%%%% Set which optimisations to perform %%%%
optimisations = {%'R_CHOI', 1, 'R(Choi)', 'high', 'SDPT3';
                 %'CAPACITY', 1, 'Capacity', 'low', 'SDPT3';
                 'CHANNEL_ROB', 1, 'R*', 'high', 'SDPT3';
                 'R_CHOI', 2, 'R_CPR', 'high', 'sed'
                 };

%%%% Amplitude damping followed by X-rotation %%%%

channel_type = 'AMP_X'; % Replace with 'X_AMP' for rotation followed by 
                        % amplitude damping
file_name = 'amp_X_results.csv';
parameter_names = {'p', 'theta'};

p = [0 0.1 ];  % List of values of the noise parameter p.

num_points = 5; % Set number of data points (different rotation angles) to
                % be calculated for each value of noise parameter p.
last_angle = pi/8; % Rotation angle will sweep from 0 to last_angle.

phase_angle = [1:num_points]*(1/num_points)*last_angle;

parameters = {p, phase_angle};




             
robustness_config.optimisations = optimisations; 
robustness_config.A_matrices = A_matrices;



amp_X_result_array = noisyChannelParamSweep(channel_type,parameter_names,...
                    parameters,robustness_config,file_name);
                


