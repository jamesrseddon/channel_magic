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
function [ choi_state ] = makeChoiState( operation )
%makeChoiState FInd Choi state for unitary or Kraus rep of channel.
%   Checks whether input is a unitary matrix or a Kraus representation of a
%   channel, then constructs the Choi state accordingly.

operation_type = checkOpType(operation);

switch operation_type
    case {'NK','NU'}
        errorStruct.message = 'Input array is not a valid CPTP channel';
        errorStruct.identifier = 'quasi:makeChoiState:notCPTP';
        error(errorStruct);
    case 'U'
        choi_state = makeJamState_U(operation);
    case 'K'
        choi_state = makeJamState_Kraus(operation);
end



