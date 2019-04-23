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
function [ kraus_operators, unitary_correction_label_ ] = ...
    pauliReset( measure_pauli_string,varargin )
%PAULIRESET Construct Kraus operators for a Pauli reset channel.
%   measure_pauli is the string identifying the Pauli to be measured.
%   After measurement, if outcome is +1, state is left as is. If outcome is
%   -1, it is corrected to a +1 eigenstate by applying a pauli that
%   anticommutes with measure_pauli. This can be specified by entering a
%   string in varargin, otherwise one will be found automatically.

% check and process inputs

if ~strcmp(class(measure_pauli_string),'char')
    error_struct.message = ['First argument should be a string '...
        'identifying the Pauli to be measured.'];
    error_struct.identifier = 'quasi:pauliReset:inputNotString';
    error(error_struct);
end

measure_pauli = pauliTensor(measure_pauli_string);
dim_measure = size(measure_pauli,1);


% Check if a unitary correction argument is supplied, and if so, check the
% input is valid and matches the dimension of the measured Pauli.

correction_specified = false;

if nargin > 1
    correction_specified = true;
    correction_pauli_string = varargin{1};
    
    if ~strcmp(class(correction_pauli_string),'char')
        error_struct.message = ['Second argument should be a string '...
            'identifying the Pauli correction to be applied.'];
        error_struct.identifier = 'quasi:pauliReset:secondInputNotString';
        error(error_struct);
    end
    
    correction_pauli = pauliTensor(correction_pauli_string);
    
    dim_correction = size(correction_pauli,1);
    
    if dim_measure ~= dim_correction
        error_struct.message = ['The dimension of the Pauli to be '...
            'measured does not match the dimension of the specified '...
            'unitary correction.'];
        error_struct.identifier = 'quasi:pauliReset:dimensionMismatch';
        error(error_struct);
    end
    
    if checkCommute(measure_pauli,correction_pauli) ~= -1
        error_struct.message = ['The correction Pauli must anticommute'...
            ' with the measurement Pauli in order to reset to +1'...
            ' eigenstate.'];
        error_struct.identifier = ['quasi:pauliReset:'...
            'inputsDoNotAnticommute'];
        error(error_struct);
    end
end

% Fix the unitary correction
if correction_specified
    unitary_correction = correction_pauli;
    unitary_correction_label = correction_pauli_string;
else
    % Find an anticommuting Pauli
    num_qubits = log2(dim_measure);
    num_paulis = 4^num_qubits;
    [pauli_array, pauli_labels] = enumeratePaulis(num_qubits);
    for kk = 1:num_paulis
        this_pauli = pauli_array(:,:,kk);
        if checkCommute(measure_pauli,this_pauli) == -1;
            unitary_correction = this_pauli;
            unitary_correction_label = pauli_labels{kk,1};
            display(['Correction applied will be '...
                unitary_correction_label '.']);
            break
        end
    end
end

% Construct the channel.
pauli_measure_kraus = pauliMeasureChannel(measure_pauli_string);

kraus_operators = conditionalChannel(pauli_measure_kraus,...
                            eye(dim_measure),unitary_correction);
                        
end


