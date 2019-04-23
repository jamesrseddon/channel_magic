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
function [ output_args ] = composeChannels(E_1,E_2)
%COMPOSECHANNELS Take two channels E_1 and E_2 and compose them: 
% E = E_2 \circ E_1. Output will be a set of Kraus operators.

op_type_1 = checkOpType(E_1);
op_type_2 = checkOpType(E_2);

if strcmp(op_type_1,'NK')|strcmp(op_type_1,'NU')
    errorStruct.message = [ 'First argument should be a complete Kraus'...
        ' representation for some channel, or a unitary matrix.'];
    errorStruct.identifier = ['quasi:composeChannels:'...
        'invalidChannel'];
    error(errorStruct);
end

if strcmp(op_type_2,'NK')|strcmp(op_type_2,'NU')
    errorStruct.message = [ 'Second argument should be a complete Kraus'...
        ' representation for some channel, or a unitary matrix.'];
    errorStruct.identifier = ['quasi:composeChannels:'...
        'invalidChannel'];
    error(errorStruct);
end

[dim_1,~,num_kraus_1] = size(E_1);
[dim_2,~,num_kraus_2] = size(E_2);

if dim_1 ~= dim_2
    errorStruct.message = [ 'Dimensions of channels do not agree.'];
    errorStruct.identifier = ['quasi:composeChannels:'...
        'mismatchDims'];
    error(errorStruct);
end



for jj = 1:num_kraus_1
    for kk = 1:num_kraus_2
        this_index = (jj-1)*num_kraus_2 + kk;
        kroutput(:,:,this_index) = E_2(:,:,kk)*E_1(:,:,jj);
    end
end

output_args = kroutput;

end

