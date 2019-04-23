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
function [ are_they_equal, difference_vector, max_difference ] =...
                    solutionChecker(solution_vector,A_matrix,...
                                mode,TARGET,varargin)
%SOLUTIONCHECKER Check solution from CVX is valid
%   Takes a vector of quasiprobabilities, calculates the corresponding
%   Pauli vector according to the supplied A_matrix, and checks that it
%   matches the target state.
%
%   The mode parameter allows a choice between checking against a target
%   state, or the Choi state for some operation.
%   MODE            DESCRIPTION
%   'state'         TARGET is a density matrix.
%   'choi'          TARGET is an operation. The Choi state will be 
%                   calculated and the solution checked against that.
%   'choi_diag'     TARGET is an operation E. The state E(|+><+|^n) will be
%                   calculated and compared with the solution.
%
%   varargin can be a tolerance for checking the solution. Otherwise
%   default tolerance will be used.

default_tolerance = 1e-10;

if nargin > 4
    check_tol = varargin{1};
else
    check_tol = default_tolerance;
end

solution_b_vector = A_matrix*solution_vector;

switch mode
    case 'state'
        target_state = TARGET;
    case 'choi'
        target_state = makeChoiState(TARGET);
    case 'choi_diag'
        target_state = makeChoiDiagState(TARGET);
end

dimension = size(target_state,1);

num_qubits = log2(dimension);

pauli_array = enumeratePaulis(num_qubits);

target_b_vector = calculateExpectationVec(target_state,pauli_array);

are_they_equal = checkEqual(target_b_vector,solution_b_vector,check_tol);

difference_vector = target_b_vector - solution_b_vector;

max_difference = max(abs(difference_vector));

end

