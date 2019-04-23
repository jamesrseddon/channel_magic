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
function tests = enumeratePaulisTest
tests = functiontests(localfunctions)
end

%%% helper functions
function boolean = checkEqualWithinTolerance(array_1,array_2,tolerance)
% Checks all elements of two arrays are equal up to some tolerance.
diff_array = abs(array_1 - array_2);
boolean_array = diff_array >= tolerance;
boolean_vector = reshape(boolean_array,1,[]);
sum_boolean = sum(boolean_vector);
boolean = sum_boolean == 0;
end

function boolean = checkCellArrayStringsEqual(cell_array_1,cell_array_2)
% Checks all elements of two cell arrays are equal - assumes cell array
% only contains strings.
boolean_array = cellfun(@strcmp,cell_array_1,cell_array_2);
num_elements = numel(boolean_array);
boolean_vector = reshape(boolean_array,1,[]);
sum_boolean = sum(boolean_vector);
boolean = sum_boolean == num_elements;
end

%%% setup and teardown
function setupOnce(testCase)

I = [1 0; 0 1];
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];

single_qubit_array(:,:,1) = I;
single_qubit_array(:,:,2) = X;
single_qubit_array(:,:,3) = Y;
single_qubit_array(:,:,4) = Z;

single_qubit_labels = {'I';'X';'Y';'Z'};

single_qubit_Z_vector = [0; 0; 0; 1];

single_qubit_Z_indices = [4];

two_qubit_array(:,:,1) = kron(I,I);
two_qubit_array(:,:,2) = kron(I,X);
two_qubit_array(:,:,3) = kron(I,Y);
two_qubit_array(:,:,4) = kron(I,Z);
two_qubit_array(:,:,5) = kron(X,I);
two_qubit_array(:,:,6) = kron(X,X);
two_qubit_array(:,:,7) = kron(X,Y);
two_qubit_array(:,:,8) = kron(X,Z);
two_qubit_array(:,:,9) = kron(Y,I);
two_qubit_array(:,:,10) = kron(Y,X);
two_qubit_array(:,:,11) = kron(Y,Y);
two_qubit_array(:,:,12) = kron(Y,Z);
two_qubit_array(:,:,13) = kron(Z,I);
two_qubit_array(:,:,14) = kron(Z,X);
two_qubit_array(:,:,15) = kron(Z,Y);
two_qubit_array(:,:,16) = kron(Z,Z);

two_qubit_labels = {'II';'IX';'IY';'IZ';'XI';'XX';'XY';'XZ';...
                        'YI';'YX';'YY';'YZ';'ZI';'ZX';'ZY';'ZZ'};

two_qubit_Z_vector = [zeros(3,1); 1; zeros(8,1); 1; zeros(2,1); 1];

two_qubit_Z_indices = [4; 13; 16];

% set up the fixture to verify output for each number of qubits n.
% Format: 'n' 'pauli array' 'pauli labels array' 'Z vector' 'Z indices'
output_fixture = {1,single_qubit_array,single_qubit_labels,...
                       single_qubit_Z_vector,single_qubit_Z_indices;
                  2,two_qubit_array,two_qubit_labels,...
                       two_qubit_Z_vector,two_qubit_Z_indices};

testCase.TestData.output_fixture = output_fixture;

% set the tolerance for checking arrays equal.
testCase.TestData.tolerance = 1e-20;
end

%%%%%%% TESTS %%%%%%%

function testCorrectArrayOutput(testCase)
% Iteratively check that the Pauli array output is as expected for
% different n, where n is number of qubits.
tolerance = testCase.TestData.tolerance;
output_fixture = testCase.TestData.output_fixture;

num_cases_to_check = size(output_fixture,1);

for kk = 1:num_cases_to_check
    num_qubits = output_fixture{kk,1};
    fixture_array = output_fixture{kk,2};
    [pauli_array,~,~,~] = enumeratePaulis(num_qubits);
    arrays_equal = checkEqualWithinTolerance(pauli_array,fixture_array,...
                                                tolerance);
    verifyTrue(testCase,arrays_equal,['The calculated Pauli '...
            'array for ' num2str(num_qubits) ...
            ' qubits was not equal to that expected.']);
end


end

% Test labels
function testCorrectLabelOutput(testCase)
% Iteratively check that the Pauli label output is as expected for
% different n, where n is number of qubits.
output_fixture = testCase.TestData.output_fixture;

num_cases_to_check = size(output_fixture,1);

for kk = 1:num_cases_to_check
    num_qubits = output_fixture{kk,1};
    label_fixture = output_fixture{kk,3};
    [~,pauli_label,~,~] = enumeratePaulis(num_qubits);
    
    
    labels_equal = checkCellArrayStringsEqual(pauli_label,label_fixture);
    
    
    verifyTrue(testCase,labels_equal,['The output list of Pauli '...
            'labels for ' num2str(num_qubits) ...
            ' qubits was not the same as expected.']);
end


end


% Test Z-Pauli functionality

function testZvectorOutput(testCase)
% Iteratively check Z vector is correct in n-qubit case for different n.
output_fixture = testCase.TestData.output_fixture;

num_cases_to_check = size(output_fixture,1);

for kk = 1:num_cases_to_check
    num_qubits = output_fixture{kk,1};
    Z_vector_fixture = output_fixture{kk,4};
    [~,~,Z_vector,~] = enumeratePaulis(num_qubits);
    
    Z_vectors_equal = isequal(Z_vector_fixture,Z_vector);
    
    verifyTrue(testCase,Z_vectors_equal,['The output Z vector'...
            ' for ' num2str(num_qubits) ...
            ' qubits was not the same as expected.']);    
end

end


function testZindicesOutput(testCase)
% Iteratively check vector of indices of Z Paulis is correct in n-qubit
% case for different n.
tolerance = testCase.TestData.tolerance;
output_fixture = testCase.TestData.output_fixture;

num_cases_to_check = size(output_fixture,1);

for kk = 1:num_cases_to_check
    num_qubits = output_fixture{kk,1};
    Z_indices_fixture = output_fixture{kk,5};
    [~,~,~,Z_indices] = enumeratePaulis(num_qubits);
    
    Z_indices_equal = checkEqualWithinTolerance(Z_indices,...
                                Z_indices_fixture,tolerance);
    verifyTrue(testCase,Z_indices_equal,['The output Z indices'...
            ' for ' num2str(num_qubits) ...
            ' qubits were not the same as expected.']);
end

end
