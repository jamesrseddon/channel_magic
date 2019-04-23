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
%%%%% make some conditional channels %%%%%%%%

% make a single_qubit Z-measure

Z_measure = pauliMeasureChannel('Z');

% stabiliser channel where Cliffords are applied based on output.

H = (1/sqrt(2))*[1 1; 1 -1]

S = [1 0; 0 1i];

HSHS = H*S*H*S;

cond_M_ZI_o1_H_o2_HSHS = conditionalChannel(Z_measure,H,HSHS);

cond_M_ZI_o1_H_o2_HSHS_CHOI = makeJamState_Kraus(cond_M_ZI_o1_H_o2_HSHS);

cond_M_ZI_o1_I_o2_H = conditionalChannel(Z_measure,eye(2),H);

cond_M_ZI_o1_I_o2_H_Choi = makeJamState_Kraus(cond_M_ZI_o1_I_o2_H);