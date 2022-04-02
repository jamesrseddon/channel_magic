function [ state_vector_set ] = ketsFromDyads( dyad_set )
%KETSFROMDYADS Extract the kets making up a set of dyads
%   Input should be of the form:
%       {ket_vec_1, bra_vec_1, dyad_1_1, ket_stabiliser_1, bra_stabiliser_1;
%       ket_vec_1, bra_vec_2, dyad_1_2, ket_stabiliser_1, bra_stabiliser_2;
%            ...}
%
%   Output will be:
%       {state_vector_1, stabiliser_list_1;
%        state_vector_2, stabiliser_list_2;
%        ...}

raw_ket_set = dyad_set(:,[1,4]);

kets_only = cell2mat(raw_ket_set(:,1)')';

[~, unique_ket_indices] = unique(kets_only,'rows');

state_vector_set = raw_ket_set(unique_ket_indices,:);



end

