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
function output = applyKraus(kraus_operators,input_operator)
%%% Apply Kraus operators to an input operator. Input should be in the
%%% form of a square matrix. Kraus operators should be input as a
%%% 3-dimensional array where the third index labels each term in the Kraus
%%% decomposition.

tolerance = 50*eps; % tolerance for checking valid Kraus decomp.

[input_dim, input_columns] = size(input_operator);

if ~(input_dim == input_columns)
    msg='Input operator must be in the form of a square matrix.';
    error(msg);
end

[kraus_dim,kraus_columns,kraus_num] = size(kraus_operators);

if ~(kraus_dim == kraus_columns)
    msg='Kraus operators must be in the form of square matrices.';
    error(msg);
end

if ~(kraus_dim == input_dim)
    msg='Kraus operator dimensions do not match input operator.';
    error(msg);
end

%%% check it's a valid CPTP map

for kk = 1:kraus_num
    K = kraus_operators(:,:,kk);
    KdaggerK(:,:,kk) = K'*K;
end

id_check = sum(KdaggerK,3);

id_diff = id_check - eye(kraus_dim);
sum_diffs = sum(sum(abs(id_diff)));

if sum_diffs > tolerance
    msg = 'K^dagger K does not sum to identity.';
    error(msg);
end

%%% Apply the channel
for jj = 1:kraus_num
    K = kraus_operators(:,:,jj);
    K_input_Kdagger(:,:,jj) = K*input_operator*K';
end

output = sum(K_input_Kdagger,3);

end

