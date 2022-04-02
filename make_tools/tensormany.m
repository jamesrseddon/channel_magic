function [tensor_product] = tensormany(matrix_list)
%TENSORMANY kron multiple matrices/vectors
%       matrix_list is a (d-by-d'-by-N)-dim array where d by d' is the
%       dimension of each matrix, and N is number of matrices.
[rows,cols,num_mat] = size(matrix_list);
tensor_product = 1;
for kk = 1:num_mat
    this_matrix = matrix_list(:,:,kk);
    tensor_product = kron(tensor_product,this_matrix);
end

