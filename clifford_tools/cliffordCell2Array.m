function [clifford_array] = cliffordCell2Array(clifford_cell)
%CLIFFORDCELL2ARRAY Convert cell array output from makeCliffordGroup to a
% 3-dimensional numerical array.

num_elements = size(clifford_cell,1);

dim = size(clifford_cell{1,2},1);

clifford_array = zeros(dim,dim,num_elements);

for kk=1:num_elements
    clifford_array(:,:,kk) = clifford_cell{kk,2};
end

end

