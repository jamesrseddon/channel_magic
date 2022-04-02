function choi_state = makePureChoiState(U)
%MAKEPURECHOISTATE Make the Choi state vector for some operator U.
%   Given operator U, make Choi state |Phi_U> = U |Phi>, where |Phi> is
%   maximally entangled Bell state:
%           |Phi> = (1/sqrt(d)) \sum_{k=1}^d |k>|k> .
%   State vector is output in standard (computational) basis.

dim = size(U,1);

if ~(dim == size(U,2))
        errorStruct.message = 'Input array is not a square matrix.';
        errorStruct.identifier = 'quasi:makePureChoiState:notSquare';
        error(errorStruct);
end

bell_state = makeBellStateDim(dim);

id = eye(dim);

U_tensor_I = kron(U,id);

choi_state = U_tensor_I*bell_state;

end

