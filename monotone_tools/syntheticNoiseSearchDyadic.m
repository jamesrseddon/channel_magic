function [results] = ...
    syntheticNoiseSearchDyadic(state_list,num_qubits,A_matrix,num_points)
%SYNTHETICNOISESEARCHDYADIC Search for states where synthetic noise confers
% an advantage by calculating dyadic negativity for different mixing 
% parameters.      

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
dyadic_neg_mat = zeros(num_states,num_points);
dyadic_neg_errors = zeros(num_states,num_points);
deltas = zeros(num_states,num_points);
status_matrix = cell(num_states,num_points);
result_summary = cell(num_states + 1,4);
result_summary(1,:) = {'Lambda+','Initial dyadic negativity',...
    'Min delta','Min delta index'};

for kk = 1:num_states
    this_state = randRho(dim);
    state_list(:,:,kk) = this_state;
    [Lambdaplus,~,~,omega_raw] = RobMagGdecomp(this_state,STABmat);
    omega = (omega_raw + omega_raw')/2;
    omega_list(:,:,kk) = omega;
    for ss = 1:num_points
        lambda = r(ss)*(Lambdaplus - 1) +1;
        noisy_state = (this_state + (lambda - 1)*omega)/lambda;
        [biglambda,~,status,cvx_slvtol] = Lambda(noisy_state,STABmat);
        status_matrix{kk,ss} = status;
        dyadic_neg_mat(kk,ss) = biglambda;
        dyadic_neg_errors(kk,ss) = cvx_slvtol;
        if ss == 1
            sampling_constant = 1/biglambda;
        end
        deltas(kk,ss) = (lambda - 1) + lambda*biglambda*sampling_constant;
    end
    [min_delta, min_delta_index] = min(deltas(kk,:));
    result_summary(kk+1,:) = {Lambdaplus,dyadic_neg_mat(kk,1),...
                                min_delta,min_delta_index};
end

results.summary = result_summary;
results.states = state_list;
results.omegas = omega_list;
results.dyadic_neg = dyadic_neg_mat;
results.dyadic_tolerance = dyadic_neg_errors;
results.status = status_matrix;
results.deltas = deltas;
end

