function [ doubled_dyad_set, doubled_A_labels, doubled_A_matrix, ...
    l1norm,reduced_dyad_set_N] = makeMinimalProductDyads( state, dyad_set, num_qubits,...
    tolerance )
%MAKEMINIMALPRODUCTDYADS  Build a minimal A matrix that makes it possible to
% find a decomposition for the 2N qubit state rho \otimes rho, by finding a
% decomposition for the N qubit state rho.

[ A_mat_N, label_cell_N, error_vector_N] = ...
    makeAmatrixDyads( dyad_set, num_qubits, tolerance);

disp(['Sum of abs errors for N-qubit construction: ' ...
    num2str(sum(abs(error_vector_N)))]);

[l1norm, ~, distrib_N, ~,~] = findRobustness(...
    state,A_mat_N,'high','SDPT3','C');

display(['N-qubit l1-norm: ' num2str(l1norm)]);

dyad_indices = find(abs(distrib_N) >= tolerance);

display(['Number of non-negligible elements: ' ...
    num2str(length(dyad_indices))]);
pause

reduced_dyad_set_N = dyad_set(dyad_indices',:);

doubled_dyad_set = makeTensorProductDyads(reduced_dyad_set_N);

[ doubled_A_matrix, doubled_A_labels, error_vector_2N] = ...
    makeAmatrixDyads( doubled_dyad_set, 2*num_qubits, tolerance);

disp(['Sum of abs errors for 2N-qubit construction: ' ...
    num2str(sum(abs(error_vector_2N)))]);

end

