function [sparse_vector,stab_indices] = sparsify(distribution,STAB_mat,N)
%SPARSIFY Sample a sparsified state vector from distribution of stabiliser
%states.
%   Input 'distribution' is the list of coefficients c_j such that:
%       |psi> = \sum_j c_j |phi_j> 
%   Where |phi_j> is the stabiliser state whose state vector is the j-th
%   column of STAB_mat. The function constructs the sparsified state:
%       |Omega> = (||c||_1/N) \sum_k^N |omega_k>
%   where each |omega_k> is chosen randomly:
%   With probability |c_j|/||c||_1, |omega_k> = c_j/|c_j| |phi_j>.
%   sparse_vector gives the state in state vector form. stab_indices lists
%   the indices of stabiliser states samples.
%
%   Requires discretesample.m, found here:
%   https://uk.mathworks.com/matlabcentral/fileexchange/
%               21912-sampling-from-a-discrete-distribution

[rows,cols] = size(STAB_mat);

abs_dist = abs(distribution);
one_norm = sum(abs_dist);
prob_dist = abs_dist/one_norm;
phase_dist = distribution./abs_dist;

stab_indices = discretesample(prob_dist,N);

sparse_vector = zeros(rows,1);

for kk = 1:N
    stab_index = stab_indices(kk);
    W_k = phase_dist(stab_index)*STAB_mat(:,stab_index);
    sparse_vector = sparse_vector + W_k;
end

sparse_vector = (one_norm/N)*sparse_vector;

end

