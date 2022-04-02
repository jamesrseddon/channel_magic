function [ state_matrix, labels ] = pauliMat2StateVec(A_matrix)
%PAULIMAT2STATEVEC Convert a matrix of Pauli vectors to corresponding
%matrix where each column is a state vector. Columns are assumed to
%correspond to pure stabiliser states.
%  

[rows, columns] = size(A_matrix);

num_qubits  = log2(rows)/2;

if mod(num_qubits,1) ~= 0
    error_struct.message = ['Number of rows does not correspond to ' ... 
                            'integer number of qubits. There should be '...
                            '4^N rows.'];
    error_struct.identifier = ['quasi:pauliMat2StateVec:'...
                                        'invalidDimension'];
    error(error_struct);
end

[pauli_array, pauli_labels] = enumeratePaulis(num_qubits);

state_matrix = zeros(2^num_qubits,columns);
labels = {};

for kk = 1:columns
    current_pauli_vec = A_matrix(:,kk);
    [current_state_density, current_state_label] = getStabiliserState(current_pauli_vec,pauli_array,...
                                            pauli_labels);
    current_state_vector = ketFromDensity(current_state_density);
    
    state_matrix(:,kk) = current_state_vector;
    labels = [labels, current_state_label];
end

end