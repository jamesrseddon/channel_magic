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
%%%% Script to make the Kraus operators for the 7-Kraus channel

psi_minus = (1/sqrt(2))*[0; 1; -1; 0];

ket00 = [1;0;0;0];

ket11 = [0;0;0;1];

ket_plus = (1/sqrt(2))*[1;1];

ket_plus2 = kron(ket_plus,ket_plus);

ket_minus = (1/sqrt(2))*[1;-1];

ket_minus2 = kron(ket_minus,ket_minus);

ket_plusY = (1/sqrt(2))*[1;1i];

ket_plusY2 = kron(ket_plusY,ket_plusY);

ket_minusY = (1/sqrt(2))*[1;-1i];

ket_minusY2 = kron(ket_minusY,ket_minusY);

K(:,:,1) = psi_minus*psi_minus';

K(:,:,2) = (1/sqrt(2))*ket00*ket00';

K(:,:,3) = (1/sqrt(2))*ket11*ket11';

K(:,:,4) = (1/sqrt(2))*ket_plus2*ket_plus2';

K(:,:,5) = (1/sqrt(2))*ket_minus2*ket_minus2';

K(:,:,6) = (1/sqrt(2))*ket_plusY2*ket_plusY2';

K(:,:,7) = (1/sqrt(2))*ket_minusY2*ket_minusY2';

for kk = 1:7
    K_dagger_K(:,:,kk) = K(:,:,kk)'*K(:,:,kk);
end

CPTP_check = sum(K_dagger_K,3);



