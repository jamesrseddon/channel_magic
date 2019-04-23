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
%%%% make single-qubit Pauli resets %%%%

I = eye(2);
X = pauliTensor('X');
Y = pauliTensor('Y');
Z = pauliTensor('Z');


X_plus_reset = conditionalChannel(pauliMeasureChannel('X'),I,Z);
X_minus_reset = conditionalChannel(pauliMeasureChannel('X'),Z,I);


Y_plus_reset = conditionalChannel(pauliMeasureChannel('Y'),I,Z);
Y_minus_reset = conditionalChannel(pauliMeasureChannel('Y'),Z,I);

Z_plus_reset = conditionalChannel(pauliMeasureChannel('Z'),I,X);
Z_minus_reset = conditionalChannel(pauliMeasureChannel('Z'),X,I);

single_qubit_pauli_resets = {};

single_qubit_pauli_resets(1,:) = {'+X', X_plus_reset};
single_qubit_pauli_resets(2,:) = {'-X', X_minus_reset};

single_qubit_pauli_resets(3,:) = {'+Y', Y_plus_reset};
single_qubit_pauli_resets(4,:) = {'-Y', Y_minus_reset};

single_qubit_pauli_resets(5,:) = {'+Z', Z_plus_reset};
single_qubit_pauli_resets(6,:) = {'-Z', Z_minus_reset};


single_qubit_pauli_resets_CHOI = {};

for kk = 1:6
    choi_state = makeJamState_Kraus(single_qubit_pauli_resets{kk,2});
    label = single_qubit_pauli_resets{kk,1};
    single_qubit_pauli_resets_CHOI(kk,:) = {label, choi_state};
end
    