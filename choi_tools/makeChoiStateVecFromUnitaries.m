function [ unitary_choi_mat ] = makeChoiStateVecFromUnitaries(unitary_cell )
%MAKECHOISTATEVECFROMUNITARIES Generate matrix of state vectors from cell
% array of unitaries.
%   Given a cell array with rows of form:
%        { '<label for k-th unitary>' U_k  }
%   where U_k is a d-by-d unitary matrix, compute the (pure) Choi state
%   vector for each U_k, and then compile the matrix where the k-th column
%   is the Choi state for U_k.

N_unitaries = size(unitary_cell,1);

dim = size(unitary_cell{1,2},1);

unitary_choi_mat = zeros(dim^2);

for kk = 1:N_unitaries;
    current_U = unitary_cell{kk,2};
    current_choi_state = makePureChoiState(current_U);
    unitary_choi_mat(:,kk) = current_choi_state;
end

end

