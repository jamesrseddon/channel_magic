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
function output = makeJamState_U(unitary)
%%% Makes the Jamiolkowski state corresponding to an input unitary matrix
U=unitary;

dim_input = size(U,1);

if ~(dim_input == size(U,2))
    msg='Input unitary not square.';
    error(msg);
end

tolerance = 1e-10;
id_check = U'*U;
id_diff = id_check - eye(dim_input);
sum_diffs = sum(sum(abs(id_diff)));

if sum_diffs > tolerance
    msg = 'Input matrix is not unitary';
    error(msg);
end

%%% Construct the extended unitary:
U_ext = kron(U,eye(dim_input));

%%% Construct the maximally entangled state ~ |00>+|11>+...+|dd>
current_ket = zeros(dim_input^2,1);

for kk = 1:dim_input
    subsystem_ket = zeros(dim_input,1);
    subsystem_ket(kk) = 1;
    this_term =  kron(subsystem_ket,subsystem_ket);
    current_ket = current_ket + this_term;
end

Phi = (1/sqrt(dim_input))*current_ket;

%%% Calculate the ket for the Jamiolkowski state:


CJket = U_ext * Phi;

%%% Calculate the density matrix

CJrho = CJket*CJket';

output = CJrho;

end

