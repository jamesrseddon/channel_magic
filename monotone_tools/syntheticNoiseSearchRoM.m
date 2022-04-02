function [results] = ...
    syntheticNoiseSearchRoM(state_list,num_qubits,A_matrix,num_points)
%SYNTHETICNOISESEARCHROM Search for states where synthetic noise confers
% an advantage by calculating RoM for different mixing parameters.
%       

%%% INPUTS AND VALIDATION %%%
[rows, cols, num_states] = size(state_list);

dim = 2^num_qubits;

if dim ~= rows || dim ~= cols
    error_struct.message = ['Dimension of input matrices does not match'...
                            ' number of qubits.'];
    error_struct.identifier = ['quasi:syntheticNoiseSearchRoM:'...
                                    'dimensionMismatch'];
    error(error_struct);
end

[STABmat,~] = pauliMat2StateVec(A_matrix);

A_cols = size(STABmat,1);

if dim ~= A_cols
    error_struct.message = ['Dimension of stabiliser state matrix does '...
                                'not match number of qubits.'];
    error_struct.identifier = ['quasi:syntheticNoiseSearchRoM:'...
                                    'dimensionMismatchAmatrix'];
    error(error_struct);
end

%%% PROCESSING %%%
r = linspace(0,1,num_points);

omega_list = zeros(dim,dim,num_states);
robustnesses = zeros(num_states,num_points);
robustness_errors = zeros(num_states,num_points);
deltas = zeros(num_states,num_points);
status_matrix = cell(num_states,num_points);
result_summary = cell(num_states + 1,4);
result_summary(1,:) = {'Lambda+','Initial robustness','Min delta',...
    'Min delta index'};

for kk = 1:num_states
    this_state = state_list(:,:,kk);
    state_list(:,:,kk) = this_state;
    [Lambdaplus,~,~,omega_raw] = RobMagGdecomp(this_state,STABmat);
    omega = (omega_raw + omega_raw')/2;
    omega_list(:,:,kk) = omega;
    for ss = 1:num_points
        lambda = r(ss)*(Lambdaplus - 1) +1;
        noisy_state = (this_state + (lambda - 1)*omega)/lambda;
        [RoMup,status,~,~,RoMlow] = ...
            findRobustness(noisy_state,A_matrix,'high','SDPT3');
        status_matrix{kk,ss} = status;
        RoM = (RoMup + RoMlow)/2;
        robustnesses(kk,ss) = RoM;
        robustness_errors(kk,ss) = (RoMup - RoMlow)/2;
        if ss == 1
            sampling_constant = 1/RoM;
        end
        deltas(kk,ss) = (lambda - 1) + lambda*RoM*sampling_constant;
    end
    [min_delta, min_delta_index] = min(deltas(kk,:));
    result_summary(kk+1,:) = {Lambdaplus,robustnesses(kk,1),...
                                min_delta,min_delta_index};
end

results.summary = result_summary;
results.states = state_list;
results.omegas = omega_list;
results.RoM = robustnesses;
results.RoMerrors = robustness_errors;
results.status = status_matrix;
results.deltas = deltas;
end

