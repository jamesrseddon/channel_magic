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
function [ dyadic_map_set ] = makeDyadicMaps(gate_set)
%MAKEDYADICMAPS Takes a gate set and enumerates corresponding dyadic maps.
%   The function takes as input a cell array describing a set of unitary
%   gates {U_j} and builds a description of the set of possible maps of the
%   form:  T(A) = U_j * A * U_k'.
%   Input cell array should be of the form:
%   gate_set = {COLUMN_OF_LABELS,COLUMN_OF_MATRICES}
%   example = {'Id', [1 0; 0 1];
%              'Z',  [1 0; 0 -1];
%              ...}
%
%   Output will be of the form:
%   dyadic_map_set = {LEFT_UNITARY_LABEL,RIGHT_UNITARY_LABEL,...
%                      LEFT_UNITARY,RIGHT_UNITARY,PAULI_LIOUVILLE_REP};
%   example = {'X','Z',...
%               [0 1;1 0],[1 0; 0 -1];
%               [0 0 i 0;
%                0 0 0 1;
%               -1 0 0 0;
%                0 1 0 0]};
%
                 

num_gates = size(gate_set,1);

map_array = {};

for jj = 1:num_gates
    left_label = gate_set{jj,1};
    left_matrix = gate_set{jj,2};
    for kk = 1:num_gates
        right_label = gate_set{kk,1};
        right_matrix = gate_set{kk,2};
        dyadic_map_array = [];
        dyadic_map_array(:,:,1) = left_matrix;
        dyadic_map_array(:,:,2) = right_matrix;
        this_PL_rep = PLdecomp(dyadic_map_array,'LR-multiply');
        this_row = {left_label,right_label,...
                    left_matrix,right_matrix,this_PL_rep};
        map_array = [map_array;
                     this_row];
    end
end
        
dyadic_map_set = map_array;

end
