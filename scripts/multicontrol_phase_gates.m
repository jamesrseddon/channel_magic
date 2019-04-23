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
%%%% Find the robustness generated for controlled-phase gates.
clear control_phase_gate affine_array A_matrix result_row cvx_result_array
clear max_up max_low results

%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%
% Affine space array and A matrix for linear system must be loaded before
% running this script. See script 'load_robustness_files.m' in root
% directory.

% flag whether to include particular types of CVX problem.
include_choi = true; % Calculate robustness of the Choi state
include_choiTP = true; % Calculate channel robustness
include_affine_search = false; % Calculate magic capacity

tic

%%% If index is j then let k = j - 1
%%% Then the controlled-phase gate for that iteration is:
%%%    diag(1,...,1,exp[i*pi/2^k])

start_index = 1;
stop_index = 12; % N.B. large k corresponds to small angle phase rotations.
                 % When normalised as R(U(pi/2^k))^k this can lead to large
                 % numerical errors.
num_qubits = 2;                 
affine_array = n_2_affine_space_array; % Replace 'n_...' as nececessary if
                                       % num_qubits is changed from 2.
A_matrix = A_mat_2; % Amend A_mat_... as necessary if num_qubits changed.


solver_precision = 'high';
solver = 'SDPT3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[plus_ket,plus_state] = makePlusTensorN(num_qubits);

results = {'index','gate'};

if include_choi
    results = [results 'Choi-R upper', 'Choi-R lower',...
    'Choi-R (normalised)', 'Choi-R lower (normalised)',...
    'normalised Choi-R estimate', 'normalised Choi-R error'];
    choi_cell = strfind(results,'normalised Choi-R estimate');
    choi_index = find(not(cellfun('isempty',choi_cell)));
end

if include_choiTP
    results = [results 'Channel-R upper', 'Channel-R lower',...
    'Channel-R (normalised)', 'Channel-R (normalised)',...
    'Channel-R estimate', 'Channel-R error'];
    channelR_cell = strfind(results,'Channel-R estimate');
    channelR_index = find(not(cellfun('isempty',channelR_cell)));
end

if include_affine_search
    results = [results, 'Affine space detail','Capacity lower',...
        'Capacity upper','Normalised capacity lower',...
        'Normalised capacity upper','Normalised capacity mean',...
        'Normalised capacity error'];
    capacity_cell = strfind(results,'Normalised capacity mean');
    capacity_index = find(not(cellfun('isempty',capacity_cell)));
end

for mm = start_index:stop_index;
    nn = mm - 1;
    control_phase_gate(:,:,mm) = nControlPhaseNqubits(nn,num_qubits);
    U = control_phase_gate(:,:,mm);
    output_state = U*plus_state*U';
    
    result_row = {nn,U};
    
    
              
    %%% Find the Choi-robustness
    if include_choi
           [choi_up,choi_status,choi_distribution,dual_var,...
               choi_low] = findRobustness(output_state,A_matrix,...
                                        solver_precision,solver);
            normalised_choi_up = choi_up^(2^nn);
            normalised_choi_low = choi_low^(2^nn);
            avg_norm_choi = mean([normalised_choi_up normalised_choi_low]);
            choi_err = (normalised_choi_up - normalised_choi_low)/2;
            result_row = [result_row choi_up, choi_low,...
                  normalised_choi_up, normalised_choi_low,...
                  avg_norm_choi, choi_err];
    end
    
    
    %%% Find TP Choi-robustness
    if include_choiTP
            [TP_choi_up,TP_choi_status,TPchoi_distribution,...
                TP_dual_var,TP_choi_low] = ...
                    findChoiCPTProbustnessDiag(U,A_matrix,...
                                        solver_precision,solver);
            normalised_TP_choi_up = TP_choi_up^(2^nn);
            normalised_TP_choi_low = TP_choi_low^(2^nn);
            avg_norm_TP_choi = ...
                mean([normalised_TP_choi_up normalised_TP_choi_low]);
            TP_choi_err = ...
                (normalised_TP_choi_up - normalised_TP_choi_low)/2;
    
            result_row = [result_row, TP_choi_up, TP_choi_low,...
                  normalised_TP_choi_up, normalised_TP_choi_low,...
                  avg_norm_TP_choi, TP_choi_err];
    end
    
    %%%% Search over affine spaces %%%%%
    if include_affine_search
        [cvx_result_array, max_low, max_up] = findAffineSpaceRobustness(...
    	                       U,affine_array,A_matrix,...
                                        solver_precision,solver);
        normalised_max_up = max_up^(2^nn);
        normalised_max_low = max_low^(2^nn);
        avg_norm = mean([normalised_max_up normalised_max_low]);
        err = (normalised_max_up - normalised_max_low)/2
        
        result_row = [ result_row, {cvx_result_array}, max_low, max_up,...
            normalised_max_low, normalised_max_up, avg_norm, err];
    end
              
    results = [results; result_row];
end

%%% plot results
    
index_vec = cell2mat(results(2:end,1));
phase_angle = pi./(2.^index_vec);

legend_cell = {};

if include_choi
    choi_R_norm = cell2mat(results(2:end,choi_index));
    choi_R_err = cell2mat(results(2:end,choi_index + 1));
    legend_cell = [legend_cell, 'Robustness of the Choi state'];
end

if include_affine_search
    max_R_norm = cell2mat(results(2:end,capacity_index));
    max_R_err = cell2mat(results(2:end,capacity_index + 1));
    legend_cell = [legend_cell, 'Magic capacity'];
end

if include_choiTP
    choiTP_norm = cell2mat(results(2:end,channelR_index));
    choiTP_err = cell2mat(results(2:end,channelR_index + 1));
    legend_cell = [legend_cell, 'Channel robustness'];
end




figure;
if include_choi
    errorbar(phase_angle,choi_R_norm,choi_R_err,'or');
    hold on;
end
if include_affine_search
    errorbar(phase_angle,max_R_norm,max_R_err,'xb');
    hold on;
end
if include_choiTP
    errorbar(phase_angle,choiTP_norm,choiTP_err,'+k');
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Phase angle, \pi/(2^k)');ylabel('R^{2^k}');
legend(legend_cell, 'Location',...
         'SouthWest');
title([num2str(num_qubits) '-qubit Clifford'...
         ' hierarchy phase-gates: diag(1,...,exp[i \pi/2^k])']);
     
%%%% Output result array %%%%
n_2_control_phase_results = results;

toc