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
function [ boolean_output ] = checkDiagonal(A)
%CHECKDIAGONAL Check whether square matrix is diagonal
%   Returns true if diagonal, false otherwise.
check_tolerance = 50*eps;

% check input

[A_rows, A_cols] = size(A);

if A_rows~=A_cols
    error_struct.message = 'Input matrix should be square.';
    error_struct.identifier = 'quasi:checkDiagonal:inputNotSquareMatrix';
    error(error_struct);
end

% pull out the diagonal of matrix A. The we check if A_diag = A.
A_diag = diag(diag(A));

is_diagonal = checkEqual(A,A_diag,check_tolerance);

boolean_output = is_diagonal;

end

