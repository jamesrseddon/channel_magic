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
function [ kraus_ops ] = conditionalChannel(measure_channel,varargin)
%CONDITIONALCHANNEL Constructs a conditional channel from building blocks.
%   Takes as input a 'measurement' channel in Kraus form. Each Kraus
%   operator is taken to represent an outcome of the measurement.
%   varargin should contain one argument for each outcome of the
%   measurement channel. These arguments should each be either a unitary or
%   a set of Kraus operators, and represents the operation that should be
%   carried out given the respective outcome.
%   EXAMPLE:
%   Let 
%   Z_measure(:,:,1) = [1 0; 0 0];
%   Z_measure(:,:,2) = [0 0; 0 1];
%
%   Z_measure is the set of Kraus operators corresponding to a measurement
%   in the single-qubit Z-basis.
%   Let U be some unitary and T be some CPTP channel, represented as a set
%   of Kraus operators.
%   
%   Then conditionalChannel(Z_measure,U,T) will output the Kraus operators
%   for the channel corresponding to the following sequence of operations:
%       1. Measure in Z-basis
%       2a. If outcome is +1, apply unitary U.
%       2b. If outcome is -1, apply CPTP map T.

% Input validation

num_channel_args = nargin-1;

m_op_type = checkOpType(measure_channel);

if ~strcmp(m_op_type,'K')
    errorStruct.message = [ 'measure_channel should be a complete Kraus'...
        ' representation for some channel, with >1 Kraus operator.'];
    errorStruct.identifier = ['quasi:conditionalChannel:'...
        'measureNotKraus'];
    error(errorStruct);
end

[rows,columns,num_outcomes] = size(measure_channel);

if num_outcomes ~= num_channel_args
    errorStruct.message = [ 'Number of channels specified does not '...
        'match number of terms in measurement channel.'];
    errorStruct.identifier = ['quasi:conditionalChannel:'...
        'wrongNumberTerms'];
    error(errorStruct);
end

for kk = 1:num_channel_args
    this_operation = varargin{kk};
    this_op_type = checkOpType(this_operation);
    if strcmp(this_op_type,'NK')||strcmp(this_op_type,'NU')
        errorStruct.message = [ 'Inputs need to be unitary matrices, or'...
           ' complete sets of Kraus operators. There is an issue with '...
           ' the ' num2str(kk) '-th input after the measurement channel.'];
        errorStruct.identifier = ['quasi:conditionalChannel:'...
            'invalidChannel'];
        error(errorStruct);
    end
    
    [op_rows,op_columns,~] = size(this_operation);
    
    if rows~=op_rows||columns~=op_columns
        errorStruct.message = [ 'Dimensions of ' num2str(kk) '-th input'...
            ' do not match dimensions of measurement channel.'];
        errorStruct.identifier = ['quasi:conditionalChannel:'...
        'dimensionMismatch'];
        error(errorStruct);
    end
end

% begin constructing new kraus operators.

kraus_ops = zeros(rows,columns,0)

for outcome = 1:num_outcomes
    % get the measurement Kraus op for this outcome.
    measure_kraus = measure_channel(:,:,outcome); 
    
    % pick out the subsequent channel for this outcome.
    conditioned_map = varargin{outcome}; 
    
    % check how many kraus operators in this channel.
    num_operators = size(conditioned_map,3); 
    for jj = 1:num_operators
        next_index = size(kraus_ops,3) + 1
        
        kraus_ops(:,:,next_index) = conditioned_map(:,:,jj)*measure_kraus;
    end
end
