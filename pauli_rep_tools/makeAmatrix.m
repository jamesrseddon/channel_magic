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
function [A_mat,labels] = makeAmatrix(state_set);
%%% Take a set of stabiliser states and output an A matrix for finding
%%% quasiprobability distribution.
%%% Expect input of the form:
%%% { 'state label 1', DENSITY_MATRIX_1;
%%%   'state label 2', DENSITY_MATRIX_2;
%%%   ...}

num_states = size(state_set,1);

first_matrix = state_set{1,2};

dim = size(first_matrix,1);

num_qubits = log2(dim);

num_paulis = 4^num_qubits;

%%% get the list of paulis
pauli_array = enumeratePaulis(num_qubits);

%%% initialise the A matrix;
A = zeros(num_paulis,num_states);
column_headers = {};

for ss = 1:num_states;
    display(['Processing state ' num2str(ss) '.']);
    current_state = state_set{ss,2};
    current_name = state_set{ss,1};
    %%% get the b vector for current state
    b_vector = calculateExpectationVec(current_state,pauli_array);
    A(:,ss) = b_vector;
    column_headers = [column_headers,current_name];
end

A_mat = A;
labels = column_headers;

