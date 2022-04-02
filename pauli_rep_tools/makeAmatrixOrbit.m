function [Amat] = makeAmatrixOrbit(alpha)
%MAKEAMATRIXORBIT Build A matrix corresponding to Clifford orbit of
% single-qubit rotation U = exp[-i*Z*alpha]
%   Builds matrix where each column is the Pauli representation of the Choi
%   state for V_{j,k} = C_j*U*C_k, where C_j and C_k are Clifford gates.

U = [exp(-1i*alpha) 0 ; 0 exp(1i*alpha)];

cliffords_U = makeCliffordGroup(1);

cliffords_U_arr = cliffordCell2Array(cliffords_U);

U_orbit = makeOrbit(U,cliffords_U_arr,10*eps,'P');

choi_array = makeChoiStateSetFromArray(U_orbit);

Amat = makeAmatrix(choi_array);
end

