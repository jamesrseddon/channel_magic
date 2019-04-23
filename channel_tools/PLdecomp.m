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
function output = PLdecomp(channel,channel_type)
%%%%%% Find the Pauli-Liouville decomposition for any channel
%%%%%% type is 'kraus','unitary', 'LR-multiply'.
%%%%%% If type is 'unitary', channel is just the unitary matrix.
%%%%%% If type is 'kraus', input should be 3-dimensional array, where third
%%%%%% index denotes each term in the Kraus decomposition.
%%%%%% eg. channel(:,:,1) would be the first Kraus operator, and so on.
%%%%%% If type is 'LR-multiply', input is a 3-dimensional array of the
%%%%%% form:
%%%%%%      channel(:,:,1) = L
%%%%%%      channel(:,:,2) = R
%%%%%% Then the map T is given by left and right-multiplying with L and R'
%%%%%% respectively, i.e.
%%%%%%      T(A) = L*A*R'  

mode = 0;

if strcmp(channel_type,'kraus')
    mode = 1;
end

if strcmp(channel_type,'unitary')
    mode = 2;
end

if strcmp(channel_type,'LR-multiply')
    mode = 3;
end

if mode == 0
    msg = ['Channel representation type "' channel_type '" not recognised'];
    error(msg);
end

switch mode
    case 1
        dim_matrix = channel(:,:,1);
    case 2
        dim_matrix = channel;
    case 3
        dim_matrix = channel(:,:,1);
end

[dim,columns] = size(dim_matrix);

if ~(dim == columns)
    msg='Input matrix not square.';
    error(msg);
end

num_qubits = log2(dim);

if ~(floor(num_qubits)==num_qubits)
    msg = ['Dimension of input matrix is ' num2str(dim) '. This is '...
            'not valid for a system of qubits.'];
    error(msg);
end

% Get the array of Paulis
[pauli_array,~] = enumeratePaulis(num_qubits);

S_rep_size = size(pauli_array,3);

S_rep = zeros(S_rep_size);

for kk = 1:S_rep_size
    %%% Calculate how channel transforms each Pauli: P -> S(P)
    this_pauli = pauli_array(:,:,kk);
    switch mode
        case 1
            %%% Apply Kraus operators
            S_pauli = applyKraus(channel,this_pauli);
        case 2
            %%% Apply unitary
            S_pauli = channel*this_pauli*channel';
        case 3
            %%% Left/right multiply with input matrices
            L = channel(:,:,1);
            R = channel(:,:,2);
            S_pauli = L*this_pauli*R';
    end
    for jj = 1:S_rep_size
        %%% calculate S_{jj,kk} = (1/d) Tr( P_jj S(P_kk))
        test_pauli = pauli_array(:,:,jj);
        S_rep(jj,kk) = (trace(test_pauli*S_pauli))/dim;
    end
end

output = S_rep;

end