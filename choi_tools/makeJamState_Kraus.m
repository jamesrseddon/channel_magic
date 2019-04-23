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
function output = makeJamState_Kraus(kraus_operators)
%%% Makes the Jamiolkowski state corresponding to an input unitary matrix
K=kraus_operators;

dim_input = size(K,1);

num_kraus = size(K,3);

if ~(dim_input == size(K,2))
    msg='Input Kraus operators not square.';
    error(msg);
end

tolerance = 1e-12;

for kk = 1:num_kraus
    K_dagger_K(:,:,kk) = K(:,:,kk)'*K(:,:,kk);
end

id_check = sum(K_dagger_K,3);
id_diff = id_check - eye(dim_input);
sum_diffs = sum(sum(abs(id_diff)));

if sum_diffs > tolerance
    msg = 'K^dagger K does not sum to identity.';
    error(msg);
end

%%% Construct the extended channel:
for kk = 1:num_kraus
    K_ext(:,:,kk) = kron(K(:,:,kk),eye(dim_input));
end

%%% Construct the maximally entangled state ~ |00>+|11>+...+|dd>
current_ket = zeros(dim_input^2,1);

for kk = 1:dim_input
    subsystem_ket = zeros(dim_input,1);
    subsystem_ket(kk) = 1;
    this_term =  kron(subsystem_ket,subsystem_ket);
    current_ket = current_ket + this_term;
end

Phi_ket = (1/sqrt(dim_input))*current_ket;

%%% Make it a density matrix
Phi = Phi_ket*Phi_ket';

%%% Calculate the density matrix for the Jamiolkowski state:
for kk=1:num_kraus
    kraus_rho(:,:,kk) = K_ext(:,:,kk) * Phi * K_ext(:,:,kk)';
end

CJrho = sum(kraus_rho,3);

output = CJrho;

end

