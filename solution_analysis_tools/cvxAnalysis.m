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
function [ output_args ] = cvxAnalysis(cvx_results,A_col_headers,tolerance)
% Gives some human-readable info on a quasiprobability distribution of a
% channel.
% A_col_headers is the set of column headers for the A matrix used in the
% linear program.
% tolerance sets the level at which we consider the contribution from a
% channel to be non-negligible.

pseudo_robustness = cvx_results{1};
quasi_distribution = cvx_results{3};    

non_zero = find(abs(quasi_distribution)>tolerance);

num_non_zero = length(non_zero);


abs_weights = 0
for kk = 1:num_non_zero;
    channel_label = A_col_headers{non_zero(kk)};
    display(['Channel: ' channel_label ', Weight = '...
        num2str(quasi_distribution(non_zero(kk)))]);
    abs_weights = abs_weights + abs(quasi_distribution(non_zero(kk)));
end

display(['Robustness relative to given set: ' num2str(pseudo_robustness)]);
display(['Number of channels in the distribution with abs(q) > '...
            num2str(tolerance) ' : ' num2str(num_non_zero)]);
display(['Sum of absolute weights of listed channels: '...
    num2str(abs_weights)]);
end

