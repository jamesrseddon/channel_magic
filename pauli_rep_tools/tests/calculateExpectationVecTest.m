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
function tests = calculateExpectationVecTest
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

%%% setup and teardown %%%

function setupOnce(testCase)
% Tolerance for array comparisons
testCase.TestData.tolerance = 10*eps;

% Single qubit matrices

I = [1 0; 0 1];
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];

single_qubit_paulis(:,:,1) = I;
single_qubit_paulis(:,:,2) = X;
single_qubit_paulis(:,:,3) = Y;
single_qubit_paulis(:,:,4) = Z;

max_mixed_one = I/2;
max_mixed_expectations = [1;0;0;0];

zero_state = [1 0; 0 0];
zero_exp = [1;0;0;1];

plus_state = 0.5*[ 1 1; 1 1];
plus_exp = [1;1;0;0];

minusY_state = 0.5*[1 1i; -1i 1];
minusY_exp = [1;0;-1;0];

T_state = 0.5*[1 (1-1i)/sqrt(2); (1+1i)/sqrt(2) 1];
T_exp = [1;1/sqrt(2);1/sqrt(2);0];

% two-qubit matrices
II = eye(4);
IX = kron(I,X);
XI = kron(X,I);
ZX = kron(Z,X);
YI = kron(Y,I);

two_qubit_paulis(:,:,1) = II;
two_qubit_paulis(:,:,2) = IX;
two_qubit_paulis(:,:,3) = XI;
two_qubit_paulis(:,:,4) = ZX;
two_qubit_paulis(:,:,5) = YI;

two_qubit_zero = kron(zero_state,zero_state);
two_qubit_zero_exp = [1;0;0;0;0];


two_qubit_zero_plus = kron(zero_state,plus_state);
two_qubit_zero_plus_exp = [1;1;0;1;0];

ZX_ent_ket = (1/sqrt(2))*(kron([1;0],(1/sqrt(2))*[1;1]) + ...
                                    kron([0;1],(1/sqrt(2))*[1;-1]));
ZX_ent_state = ZX_ent_ket*ZX_ent_ket';
ZX_ent_exp = [1;0;0;1;0];

ZX2_ent_ket = (1/sqrt(2))*(kron([0;1],(1/sqrt(2))*[1;1]) + ...
                                    kron([1;0],(1/sqrt(2))*[1;-1]));
ZX2_ent_state = ZX2_ent_ket*ZX2_ent_ket';
ZX2_ent_exp = [1;0;0;-1;0];

% HEADERS: 'Test ID','state','observables','correct output'
fixtures = {'EXP1',max_mixed_one,single_qubit_paulis,...
                                                max_mixed_expectations;
            'EXP2',zero_state,single_qubit_paulis,...
                                                zero_exp;
            'EXP3',minusY_state,single_qubit_paulis,...
                                                minusY_exp;
            'EXP4',T_state,single_qubit_paulis,...
                                                T_exp;
            'EXP5',plus_state,single_qubit_paulis,plus_exp;
            'EXP6',two_qubit_zero,two_qubit_paulis,two_qubit_zero_exp;
            'EXP7',two_qubit_zero_plus,two_qubit_paulis,...
                                                two_qubit_zero_plus_exp;
            'EXP8',ZX_ent_state,two_qubit_paulis,ZX_ent_exp;
            'EXP9',ZX2_ent_state,two_qubit_paulis,ZX2_ent_exp};
                                                    

testCase.TestData.exp_value_fixtures = fixtures;

% Fixtures for validating square Hermitian input matrix
empty_array = [];
col_vector_array = [1;0;0;0];
row_vector_array = [1 0 0 0];
rectangular_matrix = [ 1 0 0; 0 0 0; 0 0 0; 0 0 0];

iidentity = 1i*max_mixed_one;
real_not_symmetric = (1/2)*[1 -1; 1 1];

% HEADERS: 'Test ID','state','observables'
testCase.TestData.matrix_check_fixtures =...
    {'MAT1',single_qubit_paulis,single_qubit_paulis};

testCase.TestData.square_check_fixtures = ...
    {'SQMAT1',col_vector_array,single_qubit_paulis;
     'SQMAT2',row_vector_array,single_qubit_paulis;
     'SQMAT3',rectangular_matrix,single_qubit_paulis;
     'SQMAT4',rectangular_matrix',single_qubit_paulis}
 
 testCase.TestData.hermitian_check_fixtures = ...
     {'HMAT1',iidentity,single_qubit_paulis;
      'HMAT2',real_not_symmetric,single_qubit_paulis};

% Fixtures for validating dimensionality for array of observables.
too_many_dims(:,:,1,1) = I;
too_many_dims(:,:,2,1) = X;
too_many_dims(:,:,1,2) = Y;
too_many_dims(:,:,2,2) = Z;

II_rectangular = [II; 0 0 0 0];
XI_rectangular = [XI; 1 0 0 0];

rectangular_array(:,:,1) = II_rectangular;
rectangular_array(:,:,2) = XI_rectangular;

II_rectangular_2 = [II  [0; 0; 0; 0]];
XI_rectangular_2 = [XI [1; 0; 0; 0]];
rectangular_array_2(:,:,1) = II_rectangular_2;
rectangular_array_2(:,:,2) = XI_rectangular_2;

testCase.TestData.observable_ndims_check_fixtures = ...
     {'OBS_N1',zero_state,empty_array;
      'OBS_N2',zero_state,X;
      'OBS_N3',zero_state,too_many_dims};
  
testCase.TestData.observable_matchdims_check_fixtures = ...
    {'OBS_D1',two_qubit_zero,single_qubit_paulis;
     'OBS_D2',zero_state,two_qubit_paulis;
     'OBS_D3',two_qubit_zero,rectangular_array;
     'OBS_D4',two_qubit_zero,rectangular_array_2
    }
end
%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%


function testExpectationValues(testCase)
% Run through pairs of states and observable arrays, and check
% the calculated expectation values match the known values.
tol = testCase.TestData.tolerance;

fixtures = testCase.TestData.exp_value_fixtures;

num_cases_to_check = size(fixtures,1);

for kk = 1:num_cases_to_check
    state = fixtures{kk,2};
    observables = fixtures{kk,3};
    true_exp_values = fixtures{kk,4};
    test_ID = fixtures{kk,1};

    output_exp_values = calculateExpectationVec(state,observables);
    
    
    output_correct = checkEqualWithinTolerance(output_exp_values,...
                                            true_exp_values,tol);
    
    verifyTrue(testCase,output_correct,['The calculated expectation '...
        'value vector for test case ID ' test_ID ' did not match '...
        'the known values.']);
    
    if output_correct
        display(['Test case ' test_ID ' OK.']);
    end
end

end



function testInputMatrixValidation(testCase)
% Check validation for 'state' variable being a square Hermitian matrix.
fixtures = testCase.TestData.matrix_check_fixtures;

num_to_check = size(fixtures,1);

for kk = 1:num_to_check
    test_ID = fixtures{kk,1};
    state = fixtures{kk,2};
    observables = fixtures{kk,3};
    verifyError(testCase,@()calculateExpectationVec(state,observables),...
        'quasi:calculateExpectationVec:inputStateNotMatrix',...
        ['In test case ' test_ID ', input matrix validation failed as '...
                    'input had incorrect dimensionality but '...
                    'appropriate error was not thrown.']);
end

square_fixtures = testCase.TestData.square_check_fixtures;
sq_num_to_check = size(square_fixtures,1);

for kk = 1:sq_num_to_check
    test_ID = square_fixtures{kk,1};
    state = square_fixtures{kk,2};
    observables = square_fixtures{kk,3};
    verifyError(testCase,@()calculateExpectationVec(state,observables),...
        'quasi:calculateExpectationVec:inputStateNotSquare',...
        ['In test case ' test_ID ', input matrix validation failed as '...
                    'input was not square but '...
                    'appropriate error was not thrown.']);
end

h_fixtures = testCase.TestData.hermitian_check_fixtures;
h_num_to_check = size(h_fixtures,1)

for kk = 1:h_num_to_check
    test_ID = h_fixtures{kk,1};
    state = h_fixtures{kk,2};
    observables = h_fixtures{kk,3};
    verifyError(testCase,@()calculateExpectationVec(state,observables),...
        'quasi:calculateExpectationVec:inputStateNotHermitian',...
        ['In test case ' test_ID ', input matrix validation failed as '...
                    'input was not hermitian but '...
                    'appropriate error was not thrown.']);
end

end

function testSomeFunctionality(testCase)
% Check validation for 'observable' variable having correct dimensionality.
ndims_fixtures = testCase.TestData.observable_ndims_check_fixtures;

ndims_num_to_check = size(ndims_fixtures,1);

for kk = 1:ndims_num_to_check
    test_ID = ndims_fixtures{kk,1};
    state = ndims_fixtures{kk,2};
    observables = ndims_fixtures{kk,3};
    verifyError(testCase,@()calculateExpectationVec(state,observables),...
        'quasi:calculateExpectationVec:inputObservablesNot3D',...
        ['In test case ' test_ID ', validation of observables array '...
                    'failed as '...
                    'input was not a 3D array but '...
                    'appropriate error was not thrown.']);
    
end

matchdims_fixtures = testCase.TestData.observable_matchdims_check_fixtures;
matchdims_num_to_check = size(matchdims_fixtures,1);



for kk = 1:matchdims_num_to_check
    test_ID = matchdims_fixtures{kk,1};
    state = matchdims_fixtures{kk,2};
    observables = matchdims_fixtures{kk,3};
    verifyError(testCase,@()calculateExpectationVec(state,observables),...
        ['quasi:calculateExpectationVec:'...
        'inputObservablesDimensionMismatch'],...
        ['In test case ' test_ID ', validation of observables array '...
                    'failed as observable dimension '...
                    'does not match that of state, but '...
                    'appropriate error was not thrown.']);

end

end
