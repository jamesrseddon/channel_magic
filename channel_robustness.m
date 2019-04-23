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
%%%% Script to calculate channel robustness for some example channels %%%%

%%% Define some operations %%%

Tgate = [ 1 0;
               0 exp(1i*pi/4)];
           
Hadamard = (1/sqrt(2))*[1 1;
                        1 -1];
                    
% Amplitude damping noise channel with noise parameter p, Kraus operators:
%   K_1 = |0><0| + sqrt(1-p)|1><1|
%   K_2 = sqrt(p)|0><1|

amp_damp_0pt1 = amplitudeDamping(0.1);

% Amplitude damping with noise p followed by X-rotation by angle theta. 
params.p = 0.1;
params.theta = pi/8;
amp_then_Xrotate = constructChannel('AMP_X',params);

% Z measurement followed by Hadamard conditioned on "1" outcome.
Z_measure = pauliMeasureChannel('Z');
conditional_H = conditionalChannel(Z_measure,eye(2),Hadamard);

%%%% Load A matrices for linear system %%%%
load('Amat2.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_2 = A;
clear A;

load('Amat1_Clifford_Pauli_reset.mat'); % Load A matrix for CPR Choi states.

%%%% Calculate channel robustness %%%%

[R_star_Hadamard,~,distribution_Hadamard,~,~] = ...
    findChoiCPTProbustness(Hadamard,A_mat_2,'high');

[R_star_Tgate,~,distribution_Tgate,~,~] = ...
    findChoiCPTProbustness(Tgate,A_mat_2,'high');

[R_star_amp_damp,~,distribution_amp_damp,~,~] = ...
    findChoiCPTProbustness(amp_damp_0pt1,A_mat_2,'high');

[R_star_ampX,~,distribution_ampX,~,~] = ...
    findChoiCPTProbustness(amp_then_Xrotate,A_mat_2,'high');

[R_star_conditional_H,~,distribution_conditional_H,~,~] = ...
    findChoiCPTProbustness(conditional_H,A_mat_2,'high');

%%%% Calculate CPR cost %%%%

[R_CPR_Hadamard,~,distribution_CPR_Hadamard,~,~] = ...
    findChoiRobustness(Hadamard,A_mat_CP_1,'high');

[R_CPR_Tgate,~,distribution_CPR_Tgate,~,~] = ...
    findChoiRobustness(Tgate,A_mat_CP_1,'high');

[R_CPR_amp_damp,~,distribution_CPR_amp_damp,~,~] = ...
    findChoiRobustness(amp_damp_0pt1,A_mat_CP_1,'high');

[R_CPR_ampX,~,distribution_CPR_ampX,~,~] = ...
    findChoiRobustness(amp_then_Xrotate,A_mat_CP_1,'high');

[R_CPR_conditional_H,~,distribution_CPR_conditional_H,~,~] = ...
    findChoiRobustness(conditional_H,A_mat_CP_1,'high');

result_array = {'Operation','Channel robustness','R_CPR';
                'Hadamard',R_star_Hadamard,R_CPR_Hadamard;
                'T gate',R_star_Tgate,R_CPR_Tgate;
                'Amplitude damping channel',R_star_amp_damp,...
                                R_CPR_amp_damp;
                'Amplitude damped X-rotation',R_star_ampX,...
                                R_CPR_ampX;
                'Hadamard conditioned on Z measurement',...
                    R_star_conditional_H,R_CPR_conditional_H}
