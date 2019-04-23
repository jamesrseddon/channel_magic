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
function [l1norm,status,distribution,dual_var] =....
            findUrobustness(input_state,unitary,A_matrix,varargin)
% Find robustness of magic for a pure state after being acted on by
% unitary. The input state must be in column vector form.

% check inputs
dim_input = size(input_state,1);
num_columns = size(input_state,2);

if num_columns ~= 1
    msg = ['Input state must be in column vector form'];
    error(msg)
end

num_qubits = log2(dim_input);
num_paulis = 4^(num_qubits);
A_rows = size(A_matrix,1);

if num_paulis ~= A_rows
    msg = ['Expected ' num2str(num_paulis) ' Paulis but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    error(msg)
end

U=unitary;
U_dim = size(U,1);

if ~(U_dim == size(U,2))
    msg='Input unitary not square.';
    error(msg);
end

U_tolerance = 1e-12;
id_check = U'*U;
id = eye(U_dim);

if ~checkEqual(id_check,id,U_tolerance)
    msg = 'Input matrix is not unitary';
    error(msg);
end


% Calculate density matrix for the state after unitary evolution.

rho_out = unitary*input_state*input_state'*unitary';

[l1norm,status,distribution,dual_var] = findRobustness(rho_out,...
    A_matrix);
end