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
function output = deNormData(old_cell,t_column,data_column)
%%% Take a cell array of data and remove normalisation of robustness.
%%% e.g. takes input R(U(pi/2^t)^{2^t} and returns R(U(pi/2^t).

old_data = cell2mat(old_cell(2:end,data_column));

t_list = cell2mat(old_cell(2:end,t_column));

new_label = {['UNNORM: ' old_cell{1,data_column}]};

num_rows = length(old_data);

new_data = cell(num_rows,1);

for kk = 1:num_rows
    this_t = t_list(kk,1);
    new_data(kk,1) = {old_data(kk,1)^(2^(-this_t))};
end

new_column = [new_label; new_data];

output = [old_cell new_column];


