function [vector_set] = vectorizeRotationSet(matrix_set)
%VECTORIZEROTATIONSET  Vectorize set of SO(3) rotation matrices.
%   Each 3x3 matrix is converted to a 9-element vector, which becomes a
%   column of the output matrix vector_set.

num_elements = size(matrix_set,3);

vector_set = zeros(9,num_elements);

for kk = 1:num_elements
    vector_set(:,kk) = vectorizeRotation(matrix_set(:,:,kk));
end

end

