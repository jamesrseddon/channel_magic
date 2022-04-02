function [ U_array ] = Ucell2array( U_cell )
%UCELL2ARRAY Convert N-by-2 cell array where second column contains a
% unitary matrix in each entry to a d-by-d-by-N array.

if ~iscell(U_cell)
    error_struct.message = 'Input should be a cell array.';
    error_struct.identifier = 'quasi:Ucell2array:inputNotCell';
    error(error_struct)
end

test_U = U_cell{1,2};

[rows, cols] = size(test_U);

if rows ~= cols
    error_struct.message = ['Input cell array should have square '...
                        'matrices in column 2.'];
    error_struct.identifier = 'quasi:Ucell2array:invalidInput';
    error(error_struct)
end    

N = size(U_cell,1);


U_array = zeros(rows,cols,N);

for kk = 1:N
    U_k = U_cell{kk,2};
    U_array(:,:,kk) = U_k;
end

end