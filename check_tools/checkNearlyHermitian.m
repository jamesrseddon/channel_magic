function boolean = checkNearlyHermitian(matrix,tolerance)
%CHECKNEARLYHERMITIAN Does what it says on the tin.

M = matrix;
M_dagger = matrix';

boolean = checkEqual(M,M_dagger,tolerance);

end

