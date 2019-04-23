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
function [l1norm,status,distribution,dual_var] = findVectorRobustness(b_vector,A_matrix,tolerance)
%%% Get robustness for a particular state relative to a given Amatrix.
%%% State specified as a vector of Pauli expectation values.

num_paulis = size(b_vector,1);
[A_rows A_columns] = size(A_matrix);

if num_paulis ~= A_rows
    msg = ['Expected ' num2str(num_paulis) ' but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    error(msg)
end

% b_low = b_vector - tolerance;
% b_up = b_vector + tolerance;

%%% run CVX linear program
cvx_begin
    variable q(A_columns);
    dual variable y;
    minimize( norm(q,1) );
    subject to
        y :  A_matrix * q == b_vector
cvx_end

status = cvx_status;

if strcmp(status,'Solved')|strcmp(status,'Inaccurate/Solved')
    l1norm = cvx_optval;
    distribution = q;
    dual_var = y;
else
    l1norm = 0;
    distribution = 0;
    dual_var = 0;
end