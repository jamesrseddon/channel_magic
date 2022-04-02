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
function [l1norm,status,distribution,dual_var,...
    optbnd] = findRobustness(state,A_matrix,varargin)
%%% Get robustness for a particular state relative to a given Amatrix.
%%% If an argument 'high' is passed in, CVX will aim for higher precision.
%%% Second optional argument specifies the solver to be used by CVX. Pass
%%% in string 'SDPT3' to use that solver, otherwise SeDuMi will be used.
%%% Third optional argument allows the variable q to be complex.
%%%  i.e. RHO = \sum_j q_j A_j where q_j can be complex. To enable this
%%%  option, set third optional argument to 'C' or 'complex'.

dim_input = size(state,1);

num_qubits = log2(dim_input);
num_paulis = 4^(num_qubits);
[A_rows A_columns] = size(A_matrix);

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

complex_flag = false;
if num_flags > 2
    if strcmp(flags{3},'C') || strcmp(flags{3},'complex')
        complex_flag = true;
    end
end

if num_paulis ~= A_rows
    errorStruct.message = ['Expected ' num2str(num_paulis) ' but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    errorStruct.identifier = 'quasi:findRobustness:mismatchedAmatrix';
    error(errorStruct);
end

%%% calculate b vector for the input state.
pauli_array = enumeratePaulis(num_qubits);
b_vector = calculateExpectationVec(state,pauli_array,100*eps);
% Kill off imaginary part due to numerical error
% b_vector = real(b_calc_vector); 

%%% run CVX linear program.

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

if complex_flag
    solve_type = [solve_type '_complex'];
end

display(['Solve type is ' solve_type '.']);

switch solve_type
    case 'sedumi_low'
        cvx_begin
            cvx_solver sedumi
            variable q(A_columns);
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
        
    case 'sedumi_high'
        
        cvx_begin
            cvx_solver sedumi
            cvx_precision high
            variable q(A_columns);
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'SDPT_low'
        cvx_begin
            variable q(A_columns);
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'SDPT_high'
        cvx_begin
            cvx_precision high
            variable q(A_columns);
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'sedumi_low_complex'
        cvx_begin
            cvx_solver sedumi
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
        
    case 'sedumi_high_complex'
        
        cvx_begin
            cvx_solver sedumi
            cvx_precision high
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'SDPT_low_complex'
        cvx_begin
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'SDPT_high_complex'
        cvx_begin
            cvx_precision high
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
end


status = cvx_status;

if strcmp(status,'Solved')||strcmp(status,'Inaccurate/Solved')
    l1norm = cvx_optval;
    distribution = q;
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
