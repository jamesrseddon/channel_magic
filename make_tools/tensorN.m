function [tensor_product] = tensorN(A,N)
%tensorN Given matrix A, construct A^{\otimes N}

[rows,cols] = size(A);

matrix_list = zeros(rows,cols,N);

for kk = 1:N
    matrix_list(:,:,kk) = A;
end

tensor_product = tensormany(matrix_list);

end