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
function [ output ] = checkCommute(A,B)
%CHECKCOMMUTE check whether two operators commute
%   Inputs should be square matrices. 
%   If AB = BA, output is +1.
%   If AB = -BA, output is -1.
%   If they neither commute nor anticommute, output is zero.

[A_rows, A_columns] = size(A);
[B_rows, B_columns] = size(B);

% check inputs

if (A_rows ~= A_columns)|(B_rows ~= B_columns)
    error_struct.message = 'Input matrices must be square.';
    error_struct.identifier = 'quasi:checkCommute:notSquare';
    error(error_struct);
end

if A_rows ~= B_rows
    error_struct.message = 'Input matrices must be of same dimension.';
    error_struct.identifier = 'quasi:checkCommute:dimensionMismatch';
    error(error_struct);
end

% Will remain zero only in case where A,B neither commute nor anticommute.
commute_flag = 0;


zero_mat = zeros(A_rows);
check_tolerance = 50*eps;

AB = A*B;
BA = B*A;

commutator = AB - BA;

anticommutator = AB + BA;

if checkEqual(commutator,zero_mat,check_tolerance)
    commute_flag = 1;
end

if checkEqual(anticommutator,zero_mat,check_tolerance)
    commute_flag = -1;
end

output = commute_flag;

end

