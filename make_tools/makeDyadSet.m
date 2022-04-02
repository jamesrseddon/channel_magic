function [ dyad_set ] = makeDyadSet(state_set)
%MAKEDYADSET Take set of stabiliser state vectors {|PSI_j>} and construct the 
% set of dyads {|PSI_j><PSI_k|}.
%   Input should be of the form:
%       {state_vector_1, stabiliser_list_1;
%        state_vector_2, stabiliser_list_2;
%        ...}
%   Output will be the set of all combinations:
%   {ket_vec_1, bra_vec_1, dyad_1_1, ket_stabiliser_1, bra_stabiliser_1;
%    ket_vec_1, bra_vec_2, dyad_1_2, ket_stabiliser_1, bra_stabiliser_2;
%        ...}
%   

[cell_rows, cell_cols] = size(state_set);

tolerance = 100*eps

if cell_cols ~= 2 
    error_struct.message = ['Expecting cell array with two columns. '...
                                'See "help makeDyadSet" for format.'];
    error_struct.identifier = ['channel_magic:makeDyadSet:invalidInput'];
    error(error_struct)
end

dim = size(state_set{1,1},1);

for kk = 1:cell_rows
    [this_state_dim, this_state_cols] = size(state_set{kk,1});
    if this_state_cols ~= 1
        error_struct.message = ['Element in row ' num2str(kk) ' is not '...
                                    'a column vector.'];
        error_struct.identifier =['channel_magic:makeDyadSet:notVector']
        error(error_struct)
    end
    if this_state_dim ~= dim
        error_struct.message = ['Vector dimensions inconsistent in row '...
                                    num2str(kk) '.'];
        error_struct.identifier =['channel_magic:makeDyadSet:dimMismatch']
        error(error_struct)
    end
    this_state = state_set{kk,1};
    state_norm = this_state'*this_state;
    if abs(state_norm - 1) > tolerance
        error_struct.message =['State in row ' num2str(kk) ' is not '...
                                   'normalised.'];
        error_struct.identifier = ['channel_magic:makeDyadSet:'...
                                    'unnormalised'];
        error(error_struct)
    end
end                  



num_states = cell_rows;
counter = 1;
perc_complete = 0;
num_dyads = num_states^2;
dyad_set = cell(num_dyads,5);
fprintf(1,'Building %d dyads\n', num_dyads);

fprintf(1,'Dyads built: %3d%%\n', perc_complete);
for jj = 1:num_states
    ket = state_set{jj,1};
    ket_stabiliser = state_set{jj,2};
    for kk = 1:num_states
        bra = state_set{kk,1};
        bra_stabiliser = state_set{kk,2};
        dyad = ket*bra';
        new_dyad_row = {ket, bra, dyad, ket_stabiliser, bra_stabiliser};
        dyad_set(counter,:) = new_dyad_row;
        counter = counter + 1;
    end
    perc_complete = 100*counter/num_dyads;
    fprintf(1,'\b\b\b\b%3.0f%%', perc_complete);
end
fprintf(1,'\n');

end