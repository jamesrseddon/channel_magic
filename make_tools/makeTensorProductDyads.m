function [ final_dyad_set ] = makeTensorProductDyads( initial_dyad_set )
%MAKETENSORPRODUCTDYADS take a Dyad set on n qubits and construct a set of 
% tensor product dyads on 2n qubits
%   Input should be of the form:
%   {ket_vec_1, bra_vec_1, dyad_1_1, ket_stabiliser_1, bra_stabiliser_1;
%    ket_vec_1, bra_vec_2, dyad_1_2, ket_stabiliser_1, bra_stabiliser_2;
%        ...}
%   Ouput will be in same format but for 2n rather than n qubits.

num_N_dyads = size(initial_dyad_set,1);
num_2N_dyads = num_N_dyads^2;

final_dyad_set = cell(num_2N_dyads,5);

counter = 1;
perc_complete = 0;
fprintf(1,'Building %d dyads\n', num_2N_dyads);
fprintf(1,'Dyads built: %3d%%\n', perc_complete);

for jj = 1:num_N_dyads
    ket_vec_A = initial_dyad_set{jj,1};
    bra_vec_A = initial_dyad_set{jj,2};
    dyad_A = initial_dyad_set{jj,2};
    ket_stab_A = initial_dyad_set(jj,4);
    bra_stab_A = initial_dyad_set(jj,5);
    for kk = 1:num_N_dyads
    	ket_vec_B = initial_dyad_set{kk,1};
        bra_vec_B = initial_dyad_set{kk,2};
        dyad_B = initial_dyad_set{kk,2};
        ket_stab_B = initial_dyad_set(kk,4);
        bra_stab_B = initial_dyad_set(kk,5);
        stabiliser_B = initial_dyad_set(kk,4:5);
        new_ket = kron(ket_vec_A,ket_vec_B);
        new_bra = kron(bra_vec_A,bra_vec_B);
        new_dyad = kron(dyad_A,dyad_B);
        new_ket_stab = {ket_stab_A, ket_stab_B};
        new_bra_stab = {bra_stab_A, bra_stab_B};
        new_row = { new_ket, new_bra, new_dyad,...
                    new_ket_stab, new_bra_stab};
        final_dyad_set(counter,:) = new_row;
        counter = counter + 1;
    end
    perc_complete = 100*counter/num_2N_dyads;
    fprintf(1,'\b\b\b\b%3.0f%%', perc_complete);
end
fprintf(1,'\n');
end

