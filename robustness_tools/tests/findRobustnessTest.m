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
function tests = findRobustnessTest
tests = functiontests(localfunctions)
end

%%% setup and teardown
function setupOnce(testCase)
%%% Load single-qubit A matrix. Amat1.mat must be in the MATLAB path.
load('Amat1.mat');
testCase.TestData.A_mat_1 = Amat1;
%%% Load two-qubit A matrix. Amat2.mat must be in the MATLAB path.
load('Amat2.mat');
testCase.TestData.A_mat_2 = A;
%%% Load four-qubit A matrix. Amat4.mat must be in the MATLAB path.
load('Amat4.mat');
testCase.TestData.A_mat_4 = A;

%%% Generate some test states to use...
% single qubit stabiliser states:
plus = (1/2)*[1 1; 1 1];
minus = (1/2)*[1 -1; -1 1];
minusY = (1/2)*[1 1i; -1i 1]
zero = [ 1 0 ; 0 0];

testCase.TestData.single_qubit_states.zero = zero;
testCase.TestData.single_qubit_states.one = [ 0 0 ; 0 1];
testCase.TestData.single_qubit_states.plus = plus;
testCase.TestData.single_qubit_states.minus = minus;
testCase.TestData.single_qubit_states.plusY = (1/2)*[1 -1i; 1i 1];
testCase.TestData.single_qubit_states.minusY = minusY;

% single qubit magic states:
T_gate = diag([1 (1 + 1i)/sqrt(2)]);
T_state = T_gate*plus*T_gate';
testCase.TestData.single_qubit_states.T_state=T_state;

% two qubit stabiliser states:
testCase.TestData.two_qubit_states.zero_zero = kron(...
    testCase.TestData.single_qubit_states.zero,...
    testCase.TestData.single_qubit_states.zero);
testCase.TestData.two_qubit_states.plus_plus = kron(plus,plus);
testCase.TestData.two_qubit_states.minusY_minusY = kron(minusY,minusY);

bell_state_ket = (1/sqrt(2))*[1; 0; 0; 1];
bell_state = bell_state_ket*bell_state_ket';
testCase.TestData.two_qubit_states.bell_state = bell_state;


% two-qubit magic states:
T_plus = kron(T_state,plus);
testCase.TestData.two_qubit_states.T_plus_state = T_plus;
plus_T = kron(plus,T_state);
testCase.TestData.two_qubit_states.plus_T_state = plus_T;
TT_state = kron(T_state,T_state);
testCase.TestData.two_qubit_states.TT_state = TT_state;

% Choi state for the T-gate
T_tensor_I = kron(T_gate,eye(2));
T_choi_ket = T_tensor_I*bell_state;
T_choi_state = T_choi_ket*T_choi_ket';
testCase.TestData.two_qubit_states.T_choi = T_choi_state;

% four-qubit stabiliser states
[four_q_max_ent, ~] = makeMaxEntangled(4); % |Omega> maximally entangled
bell_state_product = kron(bell_state,bell_state);
plus_minusY_zero_minus = kron(plus,kron(minusY,kron(zero,minus)));

testCase.TestData.four_qubit_states.max_ent = four_q_max_ent;
testCase.TestData.four_qubit_states.plus_minusY_zero_minus = ...
                                            plus_minusY_zero_minus;
testCase.TestData.four_qubit_states.bell_product = bell_state_product;                

% four-qubit magic states
TTTT_state = kron(TT_state,TT_state);
testCase.TestData.four_qubit_states.TTTT_state = TTTT_state;

end


%%% Helper functions

%%% Test input validation
function testChecksAmatrixSize(testCase)
%%% Enter a state that has dimension not matched by the number of Paulis
%%% implied by the A matrix, and check it rejects it.

% two-qubit state and one-qubit A matrix.
two_qubit_state = testCase.TestData.two_qubit_states.zero_zero;
A_mat_1 = testCase.TestData.A_mat_1;

% one-qubit state and two-qubit A matrix.
one_qubit_state = testCase.TestData.single_qubit_states.zero;
A_mat_2 = testCase.TestData.A_mat_2;

% empty matrix
empty_thing = [];


% tests
verifyError(testCase,@()findRobustness(two_qubit_state,A_mat_1),...
        'quasi:findRobustness:mismatchedAmatrix',...
        ['Dimension check failed as state dimension does not match '...
        ' size of A matrix '...
        'but function did not throw correct exception.']);
    
verifyError(testCase,@()findRobustness(one_qubit_state,A_mat_2),...
        'quasi:findRobustness:mismatchedAmatrix',...
        ['Dimension check failed as state dimension does not match '...
        ' size of A matrix '...
        'but function did not throw correct exception.']);
    
verifyError(testCase,@()findRobustness(empty_thing,A_mat_2),...
        'quasi:findRobustness:mismatchedAmatrix',...
        ['Dimension check failed as state dimension does not match '...
        ' size of A matrix '...
        'but function did not throw correct exception.']);
    
verifyError(testCase,@()findRobustness(one_qubit_state,empty_thing),...
        'quasi:findRobustness:mismatchedAmatrix',...
        ['Dimension check failed as state dimension does not match '...
        ' size of A matrix '...
        'but function did not throw correct exception.']);
end

%%% Test stabiliser states have robustness 1
function testSingleQubitStabiliserStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_1;
error_tolerance = 1e-8; % Allowed difference between 'robustness' and 1.

% stabiliser states
state(:,:,1) = testCase.TestData.single_qubit_states.zero;
state(:,:,2) = testCase.TestData.single_qubit_states.one;
state(:,:,3) = testCase.TestData.single_qubit_states.plus;
state(:,:,4) = testCase.TestData.single_qubit_states.minus;
state(:,:,5) = testCase.TestData.single_qubit_states.plusY;
state(:,:,6) = testCase.TestData.single_qubit_states.minusY;
state_name = {'|0><0|','|1><1|','|+><+|','|-><-|','|+Y><+Y|','|-Y><-Y|'};

for kk = 1:6
    [l1norm,status,distribution,dual_var,...
        optbnd] = findRobustness(state(:,:,kk),A_matrix);
    error_size = abs(l1norm - 1);
    verifyLessThan(testCase,error_size,error_tolerance,...
        ['In test number ' num2str(kk) ', the stabiliser state '...
        state_name{kk} ' was tested, but robustness was not '...
        'calculated as equal to 1.']);
end

end

function testTwoQubitStabiliserStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_2;
error_tolerance = 1e-7; % Allowed difference between 'robustness' and 1.

% stabiliser states
state(:,:,1) = testCase.TestData.two_qubit_states.zero_zero;
state(:,:,2) = testCase.TestData.two_qubit_states.plus_plus;
state(:,:,3) = testCase.TestData.two_qubit_states.minusY_minusY;
state(:,:,4) = testCase.TestData.two_qubit_states.bell_state;
state_name = {'|00><00|','|++><++|','|-Y,-Y><-Y,-Y|','|Phi_+><Phi_+|'};

for kk = 1:4
    [l1norm,status,distribution,dual_var,...
        optbnd] = findRobustness(state(:,:,kk),A_matrix);
    error_size = abs(l1norm - 1);
    verifyLessThan(testCase,error_size,error_tolerance,...
        ['In test number ' num2str(kk) ', the stabiliser state '...
        state_name{kk} ' was tested, but robustness was not '...
        'calculated as equal to 1.']);
end

end


%%% Test T-state robustness
function testSingleQubitTStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_1;
error_tolerance = 1e-7; % Allowed difference between calculated and true R.

% State
state = testCase.TestData.single_qubit_states.T_state;
robustness = sqrt(2);

[l1norm,status,~,~,~] = findRobustness(state,A_matrix); 
error_size = abs(l1norm - robustness);
verifyLessThan(testCase,error_size,error_tolerance,...
        ['T state was tested, but the calculated robustness was not '...
        ' sqrt(2).']);
end

function testTwoQubitSingleTStateRobustness(testCase)
%%% |T>|+> should have same robustness as |T>.
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_2;
error_tolerance = 1e-7; % Allowed difference between calculated and true R.

% State
state = testCase.TestData.two_qubit_states.T_plus_state;
robustness = sqrt(2);

[l1norm,status,~,~,...
        ~] = findRobustness(state,A_matrix);
error_size = abs(l1norm - robustness);
verifyLessThan(testCase,error_size,error_tolerance,...
        ['T state was tested, but the calculated robustness was not '...
        'sqrt(2).']);
end

function testTwoQubitSingleTStateRobustnessSwitched(testCase)
%%% |+>|T> should have same robustness as |T>.
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_2;
error_tolerance = 1e-7; % Allowed difference between calculated and true R.

% State
state = testCase.TestData.two_qubit_states.plus_T_state;
robustness = sqrt(2);

[l1norm,status,distribution,~,...
        optbnd] = findRobustness(state,A_matrix);
error_size = abs(l1norm - robustness);
verifyLessThan(testCase,error_size,error_tolerance,...
        ['T state was tested, but the calculated robustness was not '...
        'sqrt(2).']);
end

function testTwoQubitTTStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_2;
error_tolerance = 1e-7; % Allowed difference between calculated and true R.

% State
state = testCase.TestData.two_qubit_states.TT_state; % |T>|T>
robustness = (1 + 3*sqrt(2))/3;

[l1norm,status,distribution,dual_var,...
        optbnd] = findRobustness(state,A_matrix);
error_size = abs(l1norm - robustness);
verifyLessThan(testCase,error_size,error_tolerance,...
        ['T tensor T state was tested, but the calculated robustness was not '...
        num2str(robustness) '.']);
end

function testTwoQubitTChoiStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_2;
error_tolerance = 1e-7; % Allowed difference between calculated and true R.

% State
state = testCase.TestData.two_qubit_states.T_choi; % (T otimes I)|Phi_+>
robustness = sqrt(2);

[l1norm,status,distribution,dual_var,...
        optbnd] = findRobustness(state,A_matrix);
error_size = abs(l1norm - robustness);
verifyLessThan(testCase,error_size,error_tolerance,...
        ['T choi state was tested, but the calculated robustness was not '...
        num2str(robustness) '.']);
end

%%% Test four-qubit stabiliser states
function testFourQubitStabiliserStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_4;
error_tolerance = 1e-6; % Allowed difference between 'robustness' and 1.

% stabiliser states
state(:,:,1) = testCase.TestData.four_qubit_states.max_ent;
state(:,:,2) = testCase.TestData.four_qubit_states.plus_minusY_zero_minus;
state(:,:,3) = testCase.TestData.four_qubit_states.bell_product;

state_name = {'|Omega><Omega|','|+,-Y,0,->','|Phi_+>|Phi_+>'};

num_to_test = size(state,3);

for kk = 1:num_to_test
    [l1norm,status,distribution,dual_var,...
        optbnd] = findRobustness(state(:,:,kk),A_matrix);
    error_size = abs(l1norm - 1);
    verifyLessThan(testCase,error_size,error_tolerance,...
        ['In test number ' num2str(kk) ', the stabiliser state '...
        state_name{kk} ' was tested, but robustness was not '...
        'calculated as equal to 1.']);
end

end


function testFourQubitTTTTStateRobustness(testCase)
%%% fixtures
% A matrix
A_matrix = testCase.TestData.A_mat_4;
error_tolerance = 1e-6; % Allowed difference between calculated and true R.

% State
state(:,:,1) = ...
    testCase.TestData.four_qubit_states.TTTT_state; % |T>|T>|T>|T>
robustness = (3 + 8 * sqrt(2))/5;
state_name = {'|T,T,T,T>'};

num_to_test = size(state,3);

for kk = 1:num_to_test
    [l1norm,status,distribution,dual_var,...
        optbnd] = findRobustness(state(:,:,kk),A_matrix);
    error_size = abs(l1norm - robustness);
    verifyLessThan(testCase,error_size,error_tolerance,...
        ['In test number ' num2str(kk) ', the stabiliser state '...
        state_name{kk} ' was tested, but robustness was not '...
        'calculated as equal to' num2str(robustness) '.']);
end

end

%%% Test other known examples
