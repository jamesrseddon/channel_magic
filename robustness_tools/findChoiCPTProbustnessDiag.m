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
    findChoiCPTProbustnessDiag(operation,A_matrix,varargin)
%%%%% find Choi-robustness, restricted to CPTP decompositions, of an input
%%%%% unitary U. This version exlploits the fact that for diagonal
%%%%% channels, R(Choi) = R(E(|+>^{\otimes n})), i.e. only n-qubit states
%%%%% need to be considered.
%%%%%
%%%%% Argument 'high' can be passed in as 3rd argument - this will pass
%%%%% through to the solver function, which will tell CVX to solve using
%%%%% the sedumi solver at high precision.


% check inputs
operation_type = checkOpType(operation);

dim_input = size(operation,1);

num_qubits = log2(dim_input);

switch operation_type
    case {'NK','NU'}
        errorStruct.message = 'Input array is not a valid CPTP channel';
        errorStruct.identifier = 'quasi:findChoiRobustnessDiag:notCPTP';
        error(errorStruct);
end


% main code

[~,plus_state] = makePlusTensorN(num_qubits);

output_state = applyKraus(operation,plus_state);

solver_precision = 'default';
solver_selection = 'default';

if nargin > 2
    solver_precision = varargin{1};
end

if nargin > 3 
    solver_selection = varargin{2};
end


[l1norm,status,distribution,dual_var,optbnd] =...
                findCPTProbustnessDiag(output_state,A_matrix,...
                    solver_precision,...
                    solver_selection);