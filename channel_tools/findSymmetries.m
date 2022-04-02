function [ symmetry_set ] = findSymmetries( A,rotation_set )
%FINDSYMMETRIES Find rotations R_j, R_k such that R_j*A*R_k = A;

num_rotations = size(rotation_set,3);
symmetry_set = [];

for jj = 1:num_rotations
    R_j = rotation_set(:,:,jj);
    for kk = 1:num_rotations
        R_k = rotation_set(:,:,kk);
        B = R_j*A*R_k;
        is_symmetry = checkEqual(A,B,1e-9);
        if is_symmetry
            symmetry_set = [symmetry_set; jj kk];
        end
    end

end

