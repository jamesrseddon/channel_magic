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
function [ A_matrix ] = makeAmatrixPLchoi(map_array,num_qubits)
%MAKEAMATRIXPLCHOI Construct A matrix from an array of PL representations of
% maps.
%   Input should be a column cell array where each cell contains the
%   Pauli-Liouville representation for a linear map.
%   Example:
%   {[1 0 0 0;
%     0 1 0 0;
%     0 0 1 0;
%     0 0 0 1];
%    [0 0 i 0;
%     0 0 0 1;
%    -1 0 0 0;
%    0 1 0 0]};
%
%    The function calculates the Pauli vector for the Choi matrix
%    corresponding to each map. In the output matrix, each column
%    corresponds to the choi matrix for a particular map.
%
%    num_qubits is input as a safety measure to ensure matrices are of the
%    correct dimension.

num_maps = length(map_array);

dim_matrix = size(map_array{1});

if ~all(dim_matrix == [4^(2*num_qubits) 4^(2*num_qubits)])
        errorStruct.message = ['Dimensions of the input Pauli-Liouville'...
                               ' matrices do not match number of qubits.'];
        errorStruct.identifier = ['quasi:makeAmatrixPLchoi:'...
                                        'dimensionMismatch'];
        error(errorStruct);
end

Phi_plus_state_PL = round(makeMaxEntangledPL(2*num_qubits));

A_matrix = [];

for kk = 1:num_maps
    map = map_array{kk};
    choi_state_PL = map*Phi_plus_state_PL;
    A_matrix = [A_matrix choi_state_PL];
end

    
end

