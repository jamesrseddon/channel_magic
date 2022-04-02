%%%%% Brute force search for equimagical decompositions %%%%%

package_folder = pwd % Should be root of repository. Set manually if 
                     % working in another directory.

load([package_folder '\robustness_files\A_comp_basis_1.mat'])
load([package_folder '\robustness_files\A_comp_basis_2.mat'])
load([package_folder '\robustness_files\Amat2_dyad.mat'])

equi_result_mat = [];


tic
y = pi/32; % initial rotation angle

initial_p = 0;
final_p = 0.2;  % noise

p_interval = 0.01;

p_vec = [initial_p:p_interval:final_p];

num_p = length(p_vec);
delta = 0.0001;

for p = p_vec
    


R_y = [ cos(2*y) -sin(2*y) 0;
        sin(2*y) cos(2*y) 0;
        0 0 1];

U_y = [exp(-1i*y) 0;
       0 exp(1i*y)];
  

plus = (1/sqrt(2))*[1;1];

psi_y = U_y * plus;


noisy_R = (1 - p)*R_y;

params.theta = y;
params.p = p;

[ kraus ] = constructChannel('Z_DEP',params);

choi_state = makeChoiState(kraus);

[initial_extent,~,~,~,~] = findPureExtent( psi_y,...
                             A_comp_basis_1,'high','SDPT3','C');

[gen_rob,gen_rob_dist] = RobMagG(choi_state,...
                                                    A_comp_basis_2);
                                          
%[cliff_dyadic,cdyad_status,cdyad_dist,~,cd_optbnd] = ...
    %findRobustness(choi_state,A_mat_1_cliffords_dyadic,'high','SDPT3','C');

[dyadic_neg,dyad_status,dyad_dist,~,dyad_optbnd] = ...
    findRobustness(choi_state,Amat2_dyad,'high','SDPT3','C');

disp(['Generalised robustness: ' num2str(gen_rob) ', '...
    'Choi state dyadic negativity: ' num2str(dyadic_neg)]);


E_0 = dyadic_neg;
E = [];
alpha = [];
Q = [];

t = 0;
prob_dist_tol = 1e-6;

continue_search = true
Amat_t=[];
while continue_search
    E_t = E_0 + delta*t;
    E = [E;E_t];
    [alpha_t, ~, ~,~ ] = angleFromExtent(E_t);
    alpha_t = real(alpha_t);
    alpha = [alpha; alpha_t]
    %Amat_cols =  makeAmatrixOrbit(alpha_t);
    %Amat_t = [Amat_t Amat_cols];
    Amat_t = makeAmatrixOrbit(alpha_t);
   
    [Q_t,status_t,dist_t,~,dist_optbnd] = ...
        findRobustness(choi_state,Amat_t,'high','SDPT3');
     p
     E_t
    pause(1)
    Q = [Q;Q_t];
    if Q_t <= 1 +prob_dist_tol
        continue_search = false;
        success = true;
    end
    if alpha_t > initial_extent
        continue_search = false;
        success = false;
    end
    t = t + 1;
end

final_T = length(Q);

if success
    disp(['Iterations: ' num2str(final_T)...
        char(10) 'Dyadic negativity: ' num2str(dyadic_neg)...
        char(10) 'Noiseless rotation extent: ' num2str(initial_extent)...
        char(10) 'Channel extent upper bound: ' num2str(E_t) ...
        char(10) 'Final Q: ' num2str(Q_t) '.']);
else

    disp(['Iterations: ' num2str(final_T)...
        char(10) 'Dyadic negativity: ' num2str(dyadic_neg)...
        char(10) 'Noiseless rotation extent: ' num2str(initial_extent)...
        char(10) 'No valid decomposition found. Final Q: ' num2str(Q_t)]);
end
toc

result_row = [y p initial_extent dyadic_neg  E_t Q_t];

equi_result_mat = [equi_result_mat; result_row];

end

result_header = {'alpha','p','noiseless extent','dyadic negativity',...
    'equimagical extent','Q'};

num_entries = length(p_vec);

equi_result_cell = [result_header;...
            mat2cell(equi_result_mat,ones(1,num_entries),ones(1,6))];
