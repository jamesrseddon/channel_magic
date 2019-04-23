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
%%%% Find the channel robustness and robustness of the Choi state for
%%%% random diagonal gates.
clear control_phase_gate affine_array A_matrix result_row cvx_result_array
clear max_up max_low results max_R_norm index_vec max_R_err choi_R_norm
clear choi_R_err random_phase_gate phase_vector
clear Rstar_err Rstar Rstar_err_idx Rstar_idx Rstar_norm Rstar_optbnd

%%% INPUTS %%%%
% Affine space array and A matrix for linear system must be loaded before
% running this script. See script 'load_robustness_files.m' in root
% directory.
affine_array = n_2_affine_space_array;  % Change for n > 2
A_matrix = A_mat_2; % Change for n > 2
num_qubits =  2; % Change for n > 2
num_points = 300; % Number of gates to randomly generate.

precision = 'high'; % Precision at which linear program terminates. Set to 
                    % 'low' for faster solution.
include_choi = true; % Calculate robustness of Choi state.
include_affine = false; % Set to true to calculate magic capacity.
% N.B. Magic capacity is much more time consuming to calculate than channel
% robustness, and we have found C = R* for the vast majority of randomly
% generated diagonal gates. Highly structured counterexamples exist, eg.
% certain multi-control phase gates - see arXiv:1901.03322.
include_Rstar = true; % Calculate channel robustness.


%%% CALCULATIONS %%%
tic
start_index = 1;
stop_index = num_points;
results = {'index'};

if include_affine
    results = [results 'C upper', 'C lower',...
    'C estimate', 'C error','affine space search detail'];
end

if include_choi
    results = [results 'Choi-R upper', 'Choi-R lower',...
    'Choi-R estimate', 'Choi-R error','Choi-R status'];
end

if include_Rstar
    results = [results 'R* lower','R* upper',...
    'R* estimate','R* error','Rstar status'];
end

results = [results,...
    'Matrix elements: <1|1>','<2|2>','<3|3>','<4|4>'];


for nn = start_index:stop_index;
    result_row = {nn};
    random_phase_gate(:,:,nn) = randomPhaseGate(num_qubits);
    this_unitary = random_phase_gate(:,:,nn);
    
    if include_affine
        [cvx_result_array, max_low, max_up] = findAffineSpaceRobustness(...
                            this_unitary,...
                            affine_array,A_matrix,precision);
        avg_max = mean([max_up max_low]);
        err = (max_up - max_low)/2;
        result_row = [result_row, max_up, max_low,...
                  avg_max, err, {cvx_result_array}];
    end
              
    %%% Find the Choi-robustness using fact that for diagonal gates,
    %%% R(E(|Omega><Omega|) = R(E(|+>^{\otimes n}))
    

    % Find robustness of Choi state
    if include_choi
        [choi_norm,choi_status,choi_distribution,dual_var,choi_optbnd] = ...
            findChoiRobustnessDiag(this_unitary,A_matrix,precision);
        choi_up = choi_norm;
        choi_low = choi_optbnd;
        avg_choi = mean([choi_up choi_low]);
        choi_err = (choi_up - choi_low)/2;
        result_row = [result_row,...
                        choi_up, choi_low,...
                        avg_choi, choi_err, choi_status];
    end
    
    % Find R*
    if include_Rstar
        [Rstar_norm,Rstar_status,Rstar_distribution,Rstar_dual_var,...
            Rstar_optbnd] = ...
            findChoiCPTProbustnessDiag(this_unitary,A_matrix,precision);
        Rstar_up = Rstar_norm;
        Rstar_low = Rstar_optbnd;
        avg_Rstar = mean([Rstar_up Rstar_low]);
        Rstar_err = (Rstar_up - Rstar_low)/2;
        
        result_row = [result_row,...
                        Rstar_up, Rstar_low,...
                        avg_Rstar, Rstar_err, Rstar_status];
    end
    
    
    phase_vector = diag(random_phase_gate(:,:,nn));
    
    result_row = [result_row,...
                  phase_vector(1),phase_vector(2),...
                  phase_vector(3),phase_vector(4)];
              
    results = [results; result_row]
end



%%% plot results


header_row = results(1,:);

if include_affine
    max_R_idx = find(strcmp([header_row],'C estimate'));
    max_R_err_idx = find(strcmp([header_row],'C error'));
    max_R = cell2mat(results(2:end,max_R_idx));
    max_R_err = cell2mat(results(2:end,max_R_err_idx));
end

if include_choi
    choi_R_idx = find(strcmp([header_row],'Choi-R estimate'));
    choi_R_err_idx = find(strcmp([header_row],'Choi-R error'));
    choi_R = cell2mat(results(2:end,choi_R_idx));
    choi_R_err = cell2mat(results(2:end,choi_R_err_idx));
end

if include_Rstar
    Rstar_idx = find(strcmp([header_row],'R* estimate'));
    Rstar_err_idx = find(strcmp([header_row],'R* error'));
    Rstar = cell2mat(results(2:end,Rstar_idx));
    Rstar_err = cell2mat(results(2:end,Rstar_err_idx));
end

diagonal_line = [ 0 5];

if include_affine && include_choi
    figure;
    plot(diagonal_line,diagonal_line,'-k');
    hold on
    errorbar(choi_R,max_R,max_R_err,'r+');
    ylabel('Magic capacity');xlabel('Robustness of the Choi state');
    title(['Capacity against robustness of Choi state'...
                ' for random diagonal gates']);
end

if include_choi && include_Rstar
    figure;
    plot(diagonal_line,diagonal_line,'-k');
    hold on
    plot(choi_R,Rstar,'r+');
    xlabel('Robustness of the Choi state');ylabel('Channel robustness');
    title(['Robustness for random diagonal gates']);
end

random_phase_results = results;

toc
