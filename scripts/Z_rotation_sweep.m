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
%%%% script for calculating measures of magic for U(theta)^{\otimes n} %%%%
%%%% where U(theta) is a single-qubit Z-rotation by angle theta %%%%

load_robustness_files;

channel_type = 'Z_tensor'
file_name = 'Z_tensor_n.csv';
parameter_names = {'n', 'theta'};
num_points = 20;
last_angle = pi/8;
phase_angle = [1:num_points]*(1/num_points)*last_angle;

N = 4; % Max number of qubits n on which to calculate R*(U^{\otimes n})
% N.B. calculating channel robustness for N=5 qubits may require
% significant amounts of memory and can take some time to run.

A_list = {A_mat_1,A_mat_2,A_mat_3,A_mat_4,A_mat_5};

results_dump = {}
for n = 1:N

    parameters = {n, phase_angle};

    A_matrices = {A_list{n}};

% 'sedumi' may be replaced with 'SDPT3' to sanity check results with a
% different solver. 
    optimisations = {'R_CHOI_DIAG', 1, 'R(Choi)', 'high', 'sedumi';
                  'CHANNEL_ROB_DIAG', 1, 'R*', 'high', 'sedumi';
                 };
             
    robustness_config.optimisations = optimisations; 
    robustness_config.A_matrices = A_matrices;


    tic
    new_results = noisyChannelParamSweep(channel_type,...
                    parameter_names,...
                    parameters,robustness_config,file_name);
    toc                


    if n == 1
        results_dump = [results_dump; new_results(1,:)];
    end
    
    results_dump = [results_dump; new_results(2:end,:)];

end

