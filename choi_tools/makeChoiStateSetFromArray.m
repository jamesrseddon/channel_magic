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
function output = makeChoiStateSetFromArray(unitary_array)
%%% Take a set of unitaries and output the corresponding set of
%%% Choi states. Expects input to be three-dimensional array where first 2
%%% dimensions label rows and columns, and third index labels the unitary
%%% element
set_size = size(unitary_array,3);
current_set = {};

for kk = 1:set_size
    disp(['Processing operation ' num2str(kk) '.']);
    jam_state = makeChoiState(unitary_array(:,:,kk));
    new_row = {num2str(kk),jam_state};
    current_set = [current_set; new_row];
end

output = current_set;