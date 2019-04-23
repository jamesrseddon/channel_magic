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
function [ out_gate_set ] = gatesetTensorIdentity( gate_set,num_id_qubits )
%GATESETTENSORIDENTITY Takes a gate set and tensors with the identity on n
% qubits.
%   Takes as input a gate set cell array of the form:
%   gate_set = {COLUMN_OF_LABELS,COLUMN_OF_MATRICES}
%   example = {'Id', [1 0; 0 1];
%              'Z',  [1 0; 0 -1];
%              ...}
%   Returns the same gate set, but with every gate U replaced with:
%       kron(U,Id), where Id is the identity for n qubits.
%

num_gates = size(gate_set,1);

dim_id = 2^num_id_qubits;

Id = eye(dim_id);

out_array = {}

for kk = 1:num_gates
    this_label = gate_set{kk,1};
    this_unitary = gate_set{kk,2};
    new_unitary = kron(this_unitary,Id);
    this_row = {this_label, new_unitary};
    out_array = [out_array; this_row];
end

out_gate_set = out_array;


end

