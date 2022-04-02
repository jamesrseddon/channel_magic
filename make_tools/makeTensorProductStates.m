function [ product_state_set] = makeTensorProductStates( initial_state_set )
%MAKETENSORPRODUCTSTATES From state set {|psi_j>} construct {|psi_j>|psi_k>}
%   Input should be:
%       {state_vector_1, stabiliser_list_1;
%        state_vector_2, stabiliser_list_2;
%        ...}
%
%   Output will be of same format but for the 2N-qubit space.

num_N_states = size(initial_state_set,1);
num_2N_states = num_N_states^2;

product_state_set = cell(num_2N_states,2);

counter = 1;
perc_complete = 0;
fprintf(1,'Building %d product states\n', num_2N_states);
fprintf(1,'States built: %3d%%\n', perc_complete);

for jj = 1:num_N_states    

    
    ket_vec_A = initial_state_set{jj,1};
    ket_stab_A = initial_state_set(jj,2);
    for kk = 1:num_N_states
    	ket_vec_B = initial_state_set{kk,1};
        ket_stab_B = initial_state_set(kk,2);
        new_ket = kron(ket_vec_A,ket_vec_B);
        new_ket_stab = {ket_stab_A, ket_stab_B};
        new_row = { new_ket, new_ket_stab};
        product_state_set(counter,:) = new_row;
        counter = counter + 1;
    end
    perc_complete = 100*(counter-1)/num_2N_states;
    fprintf(1,'\b\b\b\b%3.0f%%', perc_complete);
end
fprintf(1,'\n');
end
