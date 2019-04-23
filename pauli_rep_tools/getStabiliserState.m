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
function [state,stabiliser] = getStabiliserState(pauli_vector)
%%%% Takes vector of Pauli expectation values and calculates the stabiliser
%%%% state.
%tolerance = 1e-5;
num_paulis = size(pauli_vector,1);

num_qubits = log(num_paulis)/log(4);

dim = 2^num_qubits;

[pauli_array, pauli_labels] = enumeratePaulis(num_qubits);

stabiliser = {};

coefficient_array = reshape(repmat(pauli_vector',num_paulis,1),dim,dim,num_paulis);

signed_paulis = coefficient_array.*pauli_array;

state = (1/dim)*sum(signed_paulis,3);

for kk = 1:num_paulis
    coefficient = pauli_vector(kk);
    switch coefficient
        case 1
            current_label = pauli_labels{kk};
        case -1
            current_label = ['-' pauli_labels{kk}];
        case 0
            continue
        otherwise
            msg = 'Found a value that was not 0, +1 or -1';
            error(msg)
    end
    stabiliser = [stabiliser;current_label];
end

end
            
    

