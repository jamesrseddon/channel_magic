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
function [capacity, result_table] = findCapacity(operation,...
                                        input_matrix,A_matrix,varargin)
%%%%% Find max robustness of magic of output state after allowing operation
%%%%% to act on an input state. Maximised over a set of input states
%%%%% specified by a matrix where each column is a vector of Pauli
%%%%% expectation values for an input state.
%%%%%
%%%%% Optional arguments: 
%%%%% varargin{1}: precision; if 'high' is passed in, solves with high
%%%%%   precision.
%%%%% varargin{2}: solver; passing in 'SDPT3' will have CVX use the SDPT3
%%%%%   solver. Otherwise SeDuMi solver will be used.
%%%%% varargin{3}: filename. If this is included, the function opens the 
%%%%% specified file and writes results to this file.
%%%%%
%%%%% Input_matrix may optionally be replaced by a row vector of indices.
%%%%% In this case each index will specify a column from the A matrix used
%%%%% in the optimisation, and take this as the input state.

save_results = false;
indexed_flag = false;
precision = 'normal';
solver = 'sedumi';

if nargin > 3
    precision = varargin{1};
end

if nargin > 4
    solver = varargin{2};
end

if nargin > 5
    filename = varargin{3};
    display(['Results filename: ' filename '.']);
    save_results = true;
end

% input validation

operation_type = checkOpType(operation);

switch operation_type
    case {'NK','NU'}
        error_struct.message = 'Input array is not a valid CPTP channel';
        error_struct.identifier = 'quasi:findCapacity:notCPTP';
        error(error_struct);
end
dim_input = size(operation,1);

num_qubits = log2(dim_input);
num_paulis = 4^(2*num_qubits);
A_rows = size(A_matrix,1);


if num_paulis ~= A_rows
    error_struct.message = ['Expected ' num2str(num_paulis) ...
        ' but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    error_struct.identifier = 'quasi:findCapacity:rowCountMismatch';
    error(error_struct)
end

% Check if input is a single row - if so treat it as indices identifying
% columns of the A_matrix to use as input states. Otherwise, treat as a
% separate 
[input_rows, input_columns] = size(input_matrix);

if input_rows == 1
    indexed_flag = true;
    index_vector = input_matrix;
    input_matrix = A_matrix(:,index_vector);
    input_rows = size(A_matrix,1);
end

subgroup_size = 2^(2*num_qubits);

% Check number of rows in input vectors matches the expected number of
% Paulis for the given unitary.
if num_paulis ~= input_rows
    error_struct.message = ['Expected ' num2str(num_paulis) ...
        ' but specified input matrix '...
        'has ' num2str(input_rows) ' rows.'];
    error_struct.identifier = ['quasi:findCapacity:inputStatesRowCount'...
        'Mismatch'];
    error(error_struct)
end

result_table = {'Robustness (upper bound)', 'lower bound',...
    'difference','Input stabiliser','status','quasi distribution'};

if indexed_flag
    result_table = [result_table 'Column index'];
end

if save_results
    file_object = fopen(filename,'a');
end

% tensor with the identity
num_operators = size(operation,3);
operation_I = zeros(2*dim_input,2*dim_input,num_operators);
for mm = 1:num_operators
    operation_I(:,:,mm) = kron(operation(:,:,mm),eye(dim_input));
end

% Loop over all input states
for kk = 1:input_columns
    input_Pauli_vector = input_matrix(:,kk);
    [input_state,...
        input_stabiliser] = getStabiliserState(input_Pauli_vector);
    stabiliser_label = '{';
    for mm = 1:(subgroup_size-1)
        stabiliser_label = [stabiliser_label input_stabiliser{mm} ', '];
    end
    stabiliser_label = [stabiliser_label input_stabiliser{subgroup_size}...
                            '}'];
                        
    switch operation_type
        case 'U'
            output_state = operation_I*input_state*operation_I';
        case 'K'
            output_state = applyKraus(operation_I,input_state);
    end

    [l1norm,status,q,~,optbnd] = findRobustness(...
        output_state,A_matrix,precision,solver);
    robustness(kk) = l1norm;
    difference = l1norm - optbnd;
    result_row = {l1norm, optbnd, difference, stabiliser_label, status, q};
    
    
    
    if indexed_flag
        current_index = index_vector(kk);
        result_row = [result_row current_index];
    end
    
    result_table = [result_table;
                    result_row];
    if save_results
        add_index = '';
        if indexed_flag
            add_index = [';' num2str(current_index)];
        end
        fprintf(file_object,'%.12f;%.12f;%.12f;%s;%s%s\n',l1norm,optbnd,...
            difference,stabiliser_label,status,add_index);
    end
end

if save_results
    fclose(file_object);
end

capacity = max(robustness);

end