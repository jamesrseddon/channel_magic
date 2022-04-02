function [ extent,status,distribution,dual_var,...
    optbnd] = findUnitaryExtent(U,clifford_mat,varargin )
%FINDUNITARYEXTENT Find the extent for a unitary as defined in BBCCGH.
%   Given an n-qubit unitary U, find:
%        xi = min{(||c||_1)^2 : U = \sum_j a_j C_j}, 
%   where C_j are the n-qubit Cliffords and a_j are complex coefficients.
%   The representation used in the CVX problem setup is the pure Choi state
%   for each unitary, expressed as a state vector in the computational
%   basis. Inputs:
%       U:              The target unitary. This can either be in the
%                       standard d-by-d matrix representation, or as the 
%                       dimension d^2 state vector for the corresponding
%                       Choi state.
%                       
%       clifford_mat:   d^2-by-N matrix, where N is the number of n-qubit
%                       Cliffords. Each column is the state vector for the
%                       Choi state for a distinct Clifford. Global phase
%                       should not matter as this is taken care of by the
%                       complex coefficients {a_j}.
%
%       Option 1:       Precision. If 'high' is passed in, solve will be at
%                       higher precision.
%
%       Option 2:       Solver. If 'SDPT3' is passed in, that solver will
%                       be used. Otherwise the solver will be SeDuMi.

[U_rows, U_cols] = size(U);
cliff_rows = size(clifford_mat,1);

choi_flag = false;

solver_flag = 'default';
precision_flag = 'default';

if nargin > 2;
    precision_flag = varargin{1};
end



if nargin > 3;
    solver_flag = varargin{2};
end

if U_cols == 1
    choi_flag = true;
end

if choi_flag
    choi_dim = U_rows;
    dim = sqrt(choi_dim);
    choi_state = U;
else
    dim = U_rows;
    choi_dim = dim^2;
end

if (~choi_flag && ~(U_rows == U_cols))
    error_struct.message = ['Input array must be either a square '...
                'matrix (if using matrix rep) or a state vector'...
                '(if using Choi representation).'];
    error_struct.identifier = 'quasi:findUnitaryExtent:UisNotSquare';
    error(error_struct);
end

num_qubits = log2(dim);

if ~(mod(num_qubits,1) == 0  && num_qubits > 0)
    error_struct.message = ['Dimension of input must correspond to'...
                'an n-qubit state, i.e. dim = 2^n.'];
    error_struct.identifier = 'quasi:findUnitaryExtent:invalidDimension';
    error(error_struct);    
end

if choi_dim ~= cliff_rows
    error_struct.message = ['Dimension of input unitary does not'...
                'match that of Clifford matrix.'];
    error_struct.identifier = 'quasi:findUnitaryExtent:dimensionMismatch';
    error(error_struct);
end

if ~choi_flag
    choi_state = makePureChoiState(U);
end
    
[ extent,status,distribution,dual_var,...
    optbnd] = findPureExtent(choi_state,clifford_mat,precision_flag,...
        solver_flag);


end

