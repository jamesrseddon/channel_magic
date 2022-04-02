function [ state_set ] = Amat2STABset(A_mat)
%AMAT2STABSET Takes A matrix of Pauli vectors and outputs corresponding set
% of pure stabiliser states.
%   Output will be of the form:
%       {state_vector_1, stabiliser_list_1;
%        state_vector_2, stabiliser_list_2;
%        ...}
%       
[A_rows, A_cols] = size(A_mat);

num_qubits = log(A_rows)/log(4);

if mod(num_qubits,1) ~= 0
    errorStruct.message = ['Expecting 4^N paulis where N is number '...
                                'of qubits. N was non-integer for '...
                                'input matrix.'];
	errorStruct.identifier = ['quasi:Amat2STABset:'...
                                        'PauliCountInvalid'];
	error(errorStruct);
end

[pauli_array,pauli_labels] = enumeratePaulis(num_qubits);

result_cell = cell(A_cols, 2);

for kk = 1:A_cols
    current_A_vector = A_mat(:,kk);
    [state_rho, stabiliser_list] = getStabiliserState(current_A_vector,...
                                    pauli_array,pauli_labels);
    result_cell{kk,2} = stabiliser_list;
    state_vector = ketFromDensity(state_rho);
    result_cell{kk,1} = state_vector;
end
state_set = result_cell;
end

