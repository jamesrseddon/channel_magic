function [results,Lambda,x,sigma,omega] = ...
    syntheticNoiseSweep(rho,A_matrix,M)
%SYNTHETICNOISESWEEP Given a state of interest, compute dyadic negativities
%for a range of noisy approximations.
%   Given a state rho, an A_matrix encoding the set of stabiliser
%   states, and number of values to sweep M:
%       (1) Compute the generalised robustness and find an optimal state
%       omega for decompositions: rho = lambda*rhoprime - (lambda-1)*omega.
%       (2) Sweep through a range of mixing parameters and find the dyadic
%       negativity for the corresponding state approximation rhoprime.
%   
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

[STABmat,~] = pauliMat2StateVec(A_matrix);

stab_state_set = Amat2STABset(A_matrix);

dyad_set = makeDyadSet(stab_state_set);

[A_mat_dyadic,~,~] = makeAmatrixDyads(dyad_set,num_qubits);

[Lambda,x,sigma,omega] = RobMagGdecomp(rho,STABmat);

L = linspace(1,Lambda,M);

results = cell(M+1,5);
result_header = {'lambda', 'status','dyadic negativity','dual bound',...
    'trace check'};
results(1,:) = result_header;

for kk = 1:M
    r = L(kk);
    rhoprime = (rho + (r - 1)*omega)/r;
    trace_check = trace(rhoprime);
    [dyadic_neg,status,~,~,dual_bound] = ...
        findRobustness(rhoprime,A_mat_dyadic,'high','SDPT3','C');
    result_row = {r,status, dyadic_neg,dual_bound,trace_check};
    results(kk+1,:) = result_row;
end


end

