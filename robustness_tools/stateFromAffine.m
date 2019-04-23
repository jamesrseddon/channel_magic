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
function [output] = stateFromAffine(affine_space)
%%% Take as input an affine space specified as a binary matrix, where each
%%% row represents an element of the affine space, and output a ket vector
%%% representing the 'canonical' stabiliser state associated with that 
%%% affine space.
%%% By 'canonical' state, we mean equal superposition of all computational
%%% basis states specified by the affine space, with +1 phase on each term.
%%% Example: given the input affine space:
%%%             0 0 1 0
%%%             0 0 1 1
%%%             1 0 1 0
%%%             1 0 1 1
%%%
%%% the function will output a vector corresponding to the state:
%%%
%%% |Phi> = (|0010> + |0011> + |1010> + |1011>)/2
%%%

if ~checkBinary(affine_space)
    errorStruct.message = 'Input needs to be a binary matrix.';
    errorStruct.identifier = 'binaryTools:notBinaryMatrix';
    error(errorStruct);
end

[rows, columns] = size(affine_space);

n = columns; % number of qubits
num_elements = rows; % number of elements of the affine space.

big_N = 2^n; % Hilbert space dimension

zero_vector = zeros(big_N,1); % initialise ket with all zeroes.

state = zero_vector;

% cycle over all elements of the affine space.
for jj = 1:num_elements
    current_element = affine_space(jj,:);
    % convert binary row to a number giving the correct index in the
    % vector representation of computational basis states.
    basis_vector_index = bi2de(current_element,'left-msb') + 1;
    new_basis_ket = zero_vector;
    new_basis_ket(basis_vector_index) = 1;
    state = state + new_basis_ket; % Add the basis state to state vector. 
end

% Normalise the state.
normalised_state = (1/sqrt(num_elements))*state;

output = normalised_state;


