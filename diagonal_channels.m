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
%%%% Script to calculate magic measures for diagonal channels using affine
%%%% space techniques.

%%%% Load A matrices for linear system %%%%
load('Amat1.mat'); % Load A matrix for all single-qubit stabiliser states.
A_mat_1 = Amat1;
clear Amat1;

load('Amat2.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_2 = A;
clear A;

load('Amat3.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_3 = A;
clear A;

load('Amat4.mat'); % Load A matrix for all 2-qubit stabiliser states.
A_mat_4 = A;
clear A;

%%%% Load affine space arrays %%%%
load('affine_spaces_one_qubit.mat');
load('affine_spaces_two_qubit.mat');
load('affine_spaces_three_qubit.mat');
load('affine_spaces_four_qubit.mat');


%%%% Build some example channels %%%%

% Z-rotations

theta = pi/8;

Z_rotate_1Q = singleQubitZrot(theta);

Z_rotate_2Q = kron(Z_rotate_1Q,Z_rotate_1Q);

Z_rotate_3Q = kron(Z_rotate_2Q,Z_rotate_1Q);

Z_rotate_4Q = kron(Z_rotate_3Q,Z_rotate_1Q);


channel_list = {'single-qubit Z rotation U',Z_rotate_1Q;
                'U^{\otimes 2}',Z_rotate_2Q;
                 'U^{\otimes 3}',Z_rotate_3Q;
                  'U^{\otimes 4}',Z_rotate_4Q};

%%%% calculate magic capacity %%%              
                         
[search_result_Z1, ~,capacity_Z1] =...
    findAffineSpaceRobustness(Z_rotate_1Q,n_1_affine_space_array,A_mat_1);

[search_result_Z2, ~,capacity_Z2] =...
    findAffineSpaceRobustness(Z_rotate_2Q,n_2_affine_space_array,A_mat_2);

[search_result_Z3, ~,capacity_Z3] =...
    findAffineSpaceRobustness(Z_rotate_3Q,n_3_affine_space_array,A_mat_3);

capacity_Z4 = 0 ; % Placeholder for results table.
tic
% [search_result_Z4, ~,capacity_Z4] =...
%    findAffineSpaceRobustness(Z_rotate_4Q,n_4_affine_space_array,A_mat_4);
% Affine space search for four qubits takes a little time to run, eg around
% 12 minutes on a personal laptop. Uncomment the above lines to run this
% search.
toc

%%% Calculate channel robustness %%%

[R_star_Z1,~,distribution_Z1,~,~] =...
    findChoiCPTProbustnessDiag(Z_rotate_1Q,A_mat_1);

[R_star_Z2,~,distribution_Z2,~,~] =...
    findChoiCPTProbustnessDiag(Z_rotate_2Q,A_mat_2);

[R_star_Z3,~,distribution_Z3,~,~] =...
    findChoiCPTProbustnessDiag(Z_rotate_3Q,A_mat_3);

[R_star_Z4,~,distribution_Z4,~,~] =...
    findChoiCPTProbustnessDiag(Z_rotate_4Q,A_mat_4);



%%% RESULTS %%%


result_table_diagonal = {'Operation','Magic capacity',...
        'Channel robustness','Normalised capacity',...
        'Normalised channel robustness';
                     channel_list{1,1},capacity_Z1,R_star_Z1,...
                                         capacity_Z1,R_star_Z1;
                     channel_list{2,1},capacity_Z2,R_star_Z2,...
                                         capacity_Z2^(1/2),R_star_Z2^(1/2);
                     channel_list{3,1},capacity_Z3,R_star_Z3,...
                                         capacity_Z3^(1/3),R_star_Z3^(1/3);
                     channel_list{4,1},capacity_Z4,R_star_Z4,...
                                         capacity_Z4^(1/4),R_star_Z4^(1/4)
                                         }
