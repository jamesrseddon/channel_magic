%%%%%%%%%% Upper bound mixed channel extent for noisy-rotations %%%%%%%%%%
%%%%INPUTS%%%%
Trotate = [ 1 0;
            0 (1+1i)/sqrt(2)];
Hadamard = (1/sqrt(2))*[1 1;
                        1 -1];                  
THT = Trotate*Hadamard*Trotate;

V=THT; % Target unitary;
p_final = 0.32; % Final error probability.
p_num = 33; % Number of noise values to check.
V=THT; % Target unitary.
N = 10000; % Number of random unitary gates to generate.

%%%%LOAD PRECOMPUTED CLIFFORD DATA%%%%

package_folder = pwd % Should be root of repository. Set manually if 
                     % working in another directory.

load([package_folder '\robustness_files\A_comp_basis_2.mat'])
load([package_folder '\clifford_tools\single_qubit_Clifford_choi_states.mat'])
load([package_folder '\clifford_tools\single_qubit_Cliffords.mat'])

C_cell = clifford_1Q; % cell array of Clifford gates.
C_mat = clifford_mat; % Matrix of Clifford choi state vectors.
two_q_state_mat = A_comp_basis_2 % Matrix of two qubit stabiliser states

U_set = randomUnitarySetWithExtent( N, clifford_mat ); % Initial set of single-qubit unitaries with known 
                        % extent.


                        
tol = 0; % Tolerance for de-duping Clifford orbit.



p_vec = linspace(0,p_final,p_num);

%%%% Compile unitary set for upper bound search %%%
U_appended_set = U_set;

U_appended_set = appendCliffordOrbit(  U_set, V, C_mat,...
                             C_cell,tol);
                            
U_appended_set = appendCliffordGates(  U_appended_set, C_cell);

U_PL = PLdecomp(V, 'unitary');

xi_primal_vertices = makeChoiStateVecFromUnitaries(...
             Uarray2cell(U_appended_set.U));

extent_vec = U_appended_set.extent;

status = cell(p_num,1);
Xi = zeros(p_num,1);
witness = zeros(10,p_num);
optbnd = zeros(p_num,1);
Lambda_vec = zeros(p_num,1);
Lambda_status = cell(p_num,1);
Xi_primal = zeros(p_num,1);
Xi_primal_status = cell(p_num,1);
Xi_primal_dual = zeros(p_num,1);

for kk = 1:p_num
    
    p = p_vec(kk);
    noisy_U_PL = (1-(4*p/3))*U_PL;
    noisy_U_PL(1,1) = 1;
    
    noisy_Kraus = composeChannels(V,depolarising(p));
    noisy_Kraus_choi = makeChoiState(noisy_Kraus);
    
    [Lambda_vec(kk),~] = Lambda(noisy_Kraus_choi,...
                                                    two_q_state_mat);
    [Xi_primal(kk),Xi_primal_status{kk},~,~,...
    Xi_primal_dual(kk)] = upperBoundExtentPrimal(noisy_Kraus_choi,...
                  xi_primal_vertices,extent_vec,'high','SDPT3');

    [ SO_mat, E_vec, status{kk}, Xi(kk), witness(:,kk), optbnd(kk),...
        delta_vec, xi_vec]...
        = upperBoundChannelExtent(noisy_U_PL, U_appended_set);
    
    
end

result_matrix = [p_vec',  Xi, optbnd, Xi_primal, Xi_primal_dual,...
                    Lambda_vec];
header_row = {'p','Xi dual problem bound','dual of the dual',...
            'Xi primal problem bound','dual of primal',...
            'Dyadic negativity'};


result_cell = mat2cell(result_matrix,ones(1,p_num),...
                ones(1,size(result_matrix,2)));
            
result_tab_THT_random_10000 = result_matrix;

result_cell_THT_random_10000 = [header_row; result_cell];
