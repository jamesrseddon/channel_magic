function [ reduced_A_matrix, stabiliser_state_set ] = ...
    randomStabiliserSet(A_matrix, number_states)
% RANDOMSTABILISERSET Get a random set of stabiliser states
%   Take an A_matrix encoding the full set of stabiliser states, randomly
%   select number_states of the columns, then output the reduced A_matrix
%   and the corresponding set of stabiliser state vectors.

[A_rows, A_cols] = size(A_matrix);

if number_states > A_cols
    error_struct.message = ['Number of states specified is larger than '...
        'the number of columns in the input A matrix.'];
    error_struct.identifier = ['quasi:randomStabiliserSet:'...
        'notEnoughStates'];
    error(error_struct)
end

index_vector = randsample(A_cols,number_states);

reduced_A_matrix = A_matrix(:,index_vector);

stabiliser_state_set = Amat2STABset(reduced_A_matrix);

end

