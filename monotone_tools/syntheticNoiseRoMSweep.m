function [results,Lambda,x,sigma,omega] = ...
    syntheticNoiseRoMSweep(rho,A_matrix,M,max_param)
%SYNTHETICNOISEROMSWEEP Given a state of interest, compute robustness of 
%magic for a range of noisy approximations.
%   Given a state rho, an A_matrix encoding the set of stabiliser
%   states, and number of values to sweep M:
%       (1) Compute the generalised robustness and find an optimal state
%       omega for decompositions: rho = lambda*rhoprime - (lambda-1)*omega.
%       (2) Sweep through a range of mixing parameters and find the dyadic
%       negativity for the corresponding state approximation rhoprime.
%   max_param sets the maximum value of the mixing parameter to be used,
%   as specified by max_lambda = max_param*(Lambda - 1) + 1, i.e. as a
%   fraction of (Lambda - 1).



[A_rows,~] = size(A_matrix);
[rho_rows,~] = size(rho);
num_qubits = log2(A_rows)/2;
num_q_check = log2(rho_rows);

if num_qubits ~= num_q_check
    error_struct.message = ['Number of qubits for state does not match '...
                            'number of qubits for A matrix.'];
    error_struct.identifier = ['quasi:syntheticNoiseSweep:'...
                                    'mismatchedAmatrix'];
    error(error_struct);
end

if (max_param < 0 || max_param > 1)
    error_struct.message = ['max_param was set to ' num2str(max_param) ...
                            '. The parameter should be a number between'...
                            ' 0 and 1.'];
    error_struct.identifier = ['quasi:syntheticNoiseRoMSweep:'...
                                    'maxParamOutOfRange'];
    error(error_struct);
end                       

[STABmat,~] = pauliMat2StateVec(A_matrix);

[Lambda,x,sigma,omega] = RobMagGdecomp(rho,STABmat);

max_lambda = max_param*(Lambda - 1) +1;

L = linspace(1,max_lambda,M);

results = cell(M+1,5);
result_header = {'lambda', 'status','RoM','dual bound',...
    'trace check'};
results(1,:) = result_header;

for kk = 1:M
    r = L(kk);
    rhoprime = (rho + (r - 1)*omega)/r;
    trace_check = trace(rhoprime);
    [RoM,status,~,~,dual_bound] = ...
        findRobustness(rhoprime,A_matrix,'high','SDPT3');
    result_row = {r,status, RoM,dual_bound,trace_check};
    results(kk+1,:) = result_row;
end


end

