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
%%%%% calculate robustness for two-qubit Z-rotations
clear robustness robustness_reg rotation_angle power
clear robustness_1 robustness_reg_1 rotation_angle power
%%% inputs
start_angle = pi/8;

power = [ 0:0.25:4 ];
pauli_1 = pauliTensor('Z');
pauli_2 = pauliTensor('X');

% A_matrix_2q = A_2q_round;
A_matrix_1q = A_clifford_1q;
legend_details = {'e^{i(Z_1 + X_1)t}','e^{iZ_1t}'}

%%%%%

pauli_sum = pauli_1 + pauli_2;

%%%%%%%%%%%%%%%%%%%%%%

N = length(power);
result_list_1 ={};
result_list_2 ={};

for kk = 1:N
    rotation_angle(kk) = start_angle*10^(-power(kk));
    single_pauli_rot = singleQubitGeneralRot(rotation_angle(kk),pauli_1);
    [single_results{1},single_results{2},single_results{3},...
        single_results{4}] = unitaryQuasi(single_pauli_rot,...
        A_matrix_1q,0);
    result_list_1 = [result_list_1; single_results]
    robustness_1(kk) = single_results{1};
    robustness_reg_1(kk) = robustness_1(kk)^(1/rotation_angle(kk));
    pauli_sum_rot = singleQubitGeneralRot(rotation_angle(kk),pauli_sum)
    [cvx_results{1},cvx_results{2},cvx_results{3},...
        cvx_results{4}] = unitaryQuasi(pauli_sum_rot,A_matrix_1q,0);
    result_list_2 = [result_list_2; cvx_results];
    robustness(kk) = cvx_results{1};
    robustness_reg(kk) = robustness(kk)^(1/rotation_angle(kk));
    figure(1)
    semilogx(rotation_angle,robustness,'ob--');
    xlabel('Rotation angle');ylabel('Robustness');
    set(gca, 'xdir','reverse')
    hold on
    semilogx(rotation_angle,robustness_1,'or--');
    legend(legend_details,'Location','East');
    hold off
    drawnow
    
    figure(2)
    semilogx(rotation_angle,robustness_reg,'ob--');
    xlabel('Rotation angle');ylabel('R^{1/t}');
    set(gca, 'xdir','reverse')    
    hold on
    semilogx(rotation_angle,robustness_reg_1,'or--');
    semilogx(rotation_angle,(robustness_reg_1).^2,'k:');
    legend(legend_details,'Location','East');
    hold off
    drawnow
end
