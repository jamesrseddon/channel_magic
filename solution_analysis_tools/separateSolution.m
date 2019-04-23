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
function [positive,negative,pos_Pauli,neg_Pauli] = separateSolution(distribution,A_matrix)
%%% Takes a distribution which is the solution to the LP problem
%%% to decompose a non-stabiliser state, and separates it into
%%% positive and negative parts, i.e. in (1 + p) E_+ - p E_-, 
%%% it finds E_+ and E_-. It then calculates the corresponding
%%% Pauli vectors by multiplication with the A_matrix.

dist_pos = distribution;
dist_pos(dist_pos<0) = 0;

dist_neg = distribution;
dist_neg(dist_neg>0) = 0;

dist_neg = -dist_neg;

pos_Pauli_unnormed = A_matrix*dist_pos;
neg_Pauli_unnormed = A_matrix*dist_neg;

p = neg_Pauli_unnormed(1);

pos_Pauli = pos_Pauli_unnormed/(1+p);
neg_Pauli = neg_Pauli_unnormed/p;



positive = dist_pos/(1+p);
negative = dist_neg/p;