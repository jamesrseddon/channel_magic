function [ doubled_dyad_set, doubled_A_labels, doubled_A_matrix, ...
    l1norm,minimal_dyad_set_N] = makeSupplementedProductDyads( state,...
                dyad_set, num_qubits, supplementary_state_set, tolerance )
%MAKESUPPLEMENTEDPRODUCTDYADS  Find the minimal set of dyads that admit
% a decomposition for the 2N qubit state rho \otimes rho, by finding a
% decomposition for the N qubit state rho. Then supplement these dyads with
% additional stabiliser states. Then calculate the A matrix.

[ A_mat_N, ~, error_vector_N] = ...
    makeAmatrixDyads( dyad_set, num_qubits, tolerance);

disp(['Sum of abs errors for N-qubit construction: ' ...
    num2str(sum(abs(error_vector_N)))]);

[l1norm, ~, distrib_N, ~,~] = findRobustness(...
    state,A_mat_N,'high','SDPT3','C');

display(['N-qubit l1-norm: ' num2str(l1norm)]);

dyad_indices = find(abs(distrib_N) >= tolerance);

display(['Number of non-negligible elements: ' ...
    num2str(length(dyad_indices)) '. Press any key.']);
pause

minimal_dyad_set_N = dyad_set(dyad_indices',:);

minimal_state_set = ketsFromDyads(minimal_dyad_set_N);

product_state_set = makeTensorProductStates(minimal_state_set);

new_state_set = [product_state_set; supplementary_state_set];

display(['Total number of states: ' ...
    num2str(size(new_state_set,1)) '.']);
display(['Number of dyads to calculate: ' ...
    num2str(size(new_state_set,1)^2) '. Press any key.']);
pause
doubled_dyad_set = makeDyadSet(new_state_set);


[ doubled_A_matrix, doubled_A_labels, error_vector_2N] = ...
    makeAmatrixDyads( doubled_dyad_set, 2*num_qubits, tolerance);

disp(['Sum of abs errors for 2N-qubit construction: ' ...
    num2str(sum(abs(error_vector_2N)))]);

end

