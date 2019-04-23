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
function [l1norm,status,distribution,dual_var,optbnd] =...
                            findChoiRobustness(operation,A_matrix,varargin)
%%% Use CVX to calculate optimal Choi state decomposition for a unitary,
%%% with respect to the set of operations defined by the input A matrix.
%%% If an argument 'high' is passed in, CVX will aim for higher precision.



%%% check inputs

operation_type = checkOpType(operation);

switch operation_type
    case {'NK','NU'}
        errorStruct.message = 'Input array is not a valid CPTP channel';
        errorStruct.identifier = 'quasi:findChoiRobustness:notCPTP';
        error(errorStruct);
end
dim_input = size(operation,1);

num_qubits = log2(dim_input);
num_paulis = 4^(2*num_qubits);
A_rows = size(A_matrix,1);


if num_paulis ~= A_rows
    msg = ['Expected ' num2str(num_paulis) ' but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    error(msg)
end

%%% get Jamiolkowski (Choi) state for the unitary.
choi_state = makeChoiState(operation);

solver_precision = 'default';
solver_selection = 'default';

if nargin > 2
    solver_precision = varargin{1};
end

if nargin > 3 
    solver_selection = varargin{2};
end

%%% run CVX linear program
[l1norm,status,distribution,dual_var,optbnd] = findRobustness(choi_state,...
    A_matrix,solver_precision,solver_selection);

end

