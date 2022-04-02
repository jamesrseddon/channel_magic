function [ U_cell ] = Uarray2cell( U_array)
%UARRAY2CELL Convert d-by-d-by-N array of unitary matrices to N-by-2 cell
%array with unitaries in second column.

[rows,cols,~] = size(U_array);


if ~(ndims(U_array) == 3 && rows == cols)
    error_struct.message = ['input should be a d-by-d-by-N array, '...
                        'where d is dimension and N is number of '...
                        'unitaries in the list.'];
    error_struct.identifier = ['quasi:Uarray2cell:invalidInput'];
    error(error_struct);
end

N = size(U_array,3);

U_cell = cell(N,2);

for kk = 1:N
    U_cell{kk,1} = num2str(kk);
    U_cell{kk,2} = U_array(:,:,kk);
end

end

