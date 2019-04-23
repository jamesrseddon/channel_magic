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
function [l1norm,status,distribution,dual_var,optbnd] = ...
                    findCPTProbustnessDiag(state,A_matrix,varargin)
%%%%% find the robustness of a state, restricting allowed decompositions to
%%%%% include only mixtures of stabiliser states corresponding to
%%%%% trace-preserving maps, according to the Choi-Jamiolkowski
%%%%% isomorphism. This version of the function assumes the input is of the
%%%%% form E(|+>^{\otimes n}) where E is the channel to be decomposed. In
%%%%% this construction, the condition on trace-preserving channels is that
%%%%% in the decomp (1+p) phi_+ - p phi_-, the diagonal elements of phi_+
%%%%% and phi_- must be equal to 1/2^n. In the Pauli-Liouville picture we
%%%%% test this by checking Tr(P phi) = 0 for any pauli P comprised of only
%%%%% Zs and Is (not including the n-qubit identity).
%%%%%
%%% If an argument 'high' is passed in, CVX will aim for higher precision.


dim_input = size(state,1);

num_qubits = log2(dim_input);
num_paulis = 4^(num_qubits);
[A_rows A_columns] = size(A_matrix);

if num_paulis ~= A_rows
    error_struct.msg = ['Expected ' num2str(num_paulis)...
        ' but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    error_struct.identifier = ['quasi:findCPTProbustnessDiag:'...
                                    'dimensionMismatch'];
    error(error_struct);
end

high_precision = false;
flags = varargin(:);
num_flags = length(flags);
% Sets a flag for high_precision in solver if 'high' is the third argument.
% class(flags)
if num_flags > 0
    if strcmp(flags{1},'high')
        high_precision = true;
    end
end

SDPT_flag = false;
if num_flags > 1
    if strcmp(flags{2},'SDPT3')
        SDPT_flag = true;
    end
end

%%% calculate b vector for the input state.
[pauli_array,pauli_labels,Z_flags,Z_indices] = enumeratePaulis(num_qubits);
b_calc_vector = calculateExpectationVec(state,pauli_array);
b_vector = real(b_calc_vector);


%%% Set up the matrices to be used in CVX.
A_prime = sparse([A_matrix -A_matrix]);
A_plus = sparse([A_matrix sparse(A_rows,A_columns)]);

%%% Construct the M matrices for Paulis containing only Zs.
M_rows = 2^num_qubits - 1;
M_columns = num_paulis;
M = sparse(M_rows,M_columns);
for kk = 1:M_rows
    M(kk,Z_indices(kk))=1;
end

M_A_product =  M * A_plus;

%p_zeroes = zeros(2*A_columns,1);
%TP_zeroes = zeros(M_rows);

% Set the solver precision and selection

solve_type = 'sedumi_low';

if high_precision && ~SDPT_flag
    solve_type = 'sedumi_high';
end

if ~high_precision && SDPT_flag
    solve_type = 'SDPT_low';
end

if high_precision && SDPT_flag
    solve_type = 'SDPT_high';
end

display(['Solve type is ' solve_type '.']);

%%% run CVX problem
switch solve_type
    case 'sedumi_low'
        cvx_begin
            cvx_solver sedumi;
            variable p(2*A_columns) nonnegative;
            dual variable y;
            minimize( norm(p,1) );
            subject to
                y :  A_prime * p == b_vector
                M_A_product * p == 0
        cvx_end
    case 'sedumi_high'
        cvx_begin
            cvx_solver sedumi;
            cvx_precision high;
            variable p(2*A_columns) nonnegative;
            dual variable y;
            minimize( norm(p,1) );
            subject to
                y :  A_prime * p == b_vector
                M_A_product * p == 0
        cvx_end   
    case 'SDPT_low'
        cvx_begin
            cvx_solver sdpt3;
            variable p(2*A_columns) nonnegative;
            dual variable y;
            minimize( norm(p,1) );
            subject to
                y :  A_prime * p == b_vector
                M_A_product * p == 0
        cvx_end  
    case 'SDPT_high'
        cvx_begin
            cvx_solver sdpt3;
            cvx_precision high;
            variable p(2*A_columns) nonnegative;
            dual variable y;
            minimize( norm(p,1) );
            subject to
                y :  A_prime * p == b_vector
                M_A_product * p == 0
        cvx_end  
end

status = cvx_status;

if strcmp(status,'Solved')|strcmp(status,'Inaccurate/Solved')
    l1norm = cvx_optval;
    distribution = p;
    dual_var = y;
    if exist('cvx_optbnd');
        optbnd = cvx_optbnd;
    elseif exist('cvx_bound')
        optbnd = cvx_bound;
    end
else
    l1norm = 0;
    distribution = 0;
    dual_var = 0;
    optbnd = 0;
end
