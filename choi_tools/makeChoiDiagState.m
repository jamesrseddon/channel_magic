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
function [ choi_diag_state ] = makeChoiDiagState( operation )
%MAKECHOIDIAGSTATE Makes n-qubit state isomorphic to n-qubit diagonal
%channel.
%   Given an n-qubit diagonal channel E, constructs the n-qubit state
%   E(|+><+|^{\otimes n}), as a density matrix.

operation_type = checkOpType(operation);

switch operation_type
    case {'NK','NU'}
        errorStruct.message = 'Input array is not a valid CPTP channel';
        errorStruct.identifier = 'quasi:makeChoiState:notCPTP';
        error(errorStruct);
    case {'U','K'}
        dimension = size(operation,1);
        num_qubits = log2(dimension);
        [~,plus_state] = makePlusTensorN(num_qubits);
        choi_diag_state = applyKraus(operation,plus_state);
end


end

