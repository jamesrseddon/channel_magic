% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
%  Copyright © 2019 James R. Seddon
%  This file is part of project: Channel Magic
%
function [ valid, l1norm, difference ] = checkSolution(distribution,...
    A_matrix,state,varargin)
%CHECKSOLUTION check quasibprobability distribution against original state.
%   [valid,norm,D] = CHECKSOLUTION(q,A,rho,...) Takes a distribution, in
%   the form of a vector of coefficients q, corresponding to columns in the
%   matrix of Pauli vectors A. A*q produces the Pauli vector corresponding
%   to the distribution. This is compared to the Pauli decomposition of the
%   state rho, which is input as a density matrix. Output valid = true if
%   the vectors are equal within default tolerance 50*eps, false otherwise.
%   the output norm is the l1-norm of the distribution. Output D is the
%   difference between A*q and the correct Pauli vector.
%   
%   [valid,norm,D] = CHECKSOLUTION(q,A,rho,tol) As above, but the
%   tolerance for the comparison is set by the argument tol.

default_tol = 50*eps;

% input validation

if nargin > 3
    tol = varargin{1};
else
    tol = default_tol;
end

if ~isnumeric(tol)
    error_struct.message = 'Input tolerance should be numeric.';
    error_struct.identifier = 'quasi:checkSolution:toleranceNotNumeric';
    error(error_struct)
end

[dist_rows,dist_cols] = size(distribution);

if dist_cols ~= 1
    error_struct.message = ['Input distribution should be in the form '... 
        'of a column vector.'];
    error_struct.identifier = 'quasi:checkSolution:distributionNotVector';
    error(error_struct)    
end

[A_rows,A_cols] = size(A_matrix);

if A_cols ~= dist_rows 
    error_struct.message = ['Number of columns of the A matrix is '...
        num2str(A_cols) ', but the number of rows of the distribution '...
        'vector is ' num2str(dist_rows) '. They should be the same.'];
    error_struct.identifier = 'quasi:checkSolution:dimensionMismatch';
    error(error_struct)    
end

[state_rows,state_cols] = size(state);

if state_rows ~= state_cols
    error_struct.message = ['State should be input as density matrix, '...
        'therefore should be square matrix.'];
    error_struct.identifier = 'quasi:checkSolution:inputStateNotSquare';
    error(error_struct)    
end

if ~checkNearlyHermitian(state,50*eps)
    error_struct.message = ['State should be input as density matrix, '...
        'therefore should be Hermitian.'];
    error_struct.identifier = 'quasi:checkSolution:inputStateNotHermitian';
    error(error_struct)    
end

num_qubits = log2(state_rows);

if floor(num_qubits) ~= num_qubits
    error_struct.message = ['Not an n-qubit state. The density matrix '...
        'has dimension ' num2str(state_rows) ', which is not a power '...
        'of 2.'];
    error_struct.identifier = 'quasi:checkSolution:inputStateNotnQubit';
    error(error_struct)    
end

num_paulis = 4^(num_qubits);

if num_paulis ~= A_rows
    error_struct.message = ['For a ' num2str(num_qubits) '-qubit '...
        'system, there are ' num2str(num_paulis) ' Pauli matrices, '...
        'but the input A matrix has ' num2str(A_rows) ' rows.'];
    error_struct.identifier = 'quasi:checkSolution:wrongNumArows';
    error(error_struct)
end

%%% calculate the correct Pauli vector for input state.
pauli_array = enumeratePaulis(num_qubits);
true_pauli_vector = calculateExpectationVec(state,pauli_array);

%%% calulate Pauli vector corresponding to distribution to be checked.
dist_pauli_vector = A_matrix*distribution;

%%% compare vectors and calculate outputs.
valid = checkEqual(true_pauli_vector,dist_pauli_vector,tol);

l1norm = sum(abs(distribution));

difference = true_pauli_vector - dist_pauli_vector;

end


