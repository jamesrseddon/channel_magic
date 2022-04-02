function [ unitary_set ] = randomUnitarySetWithExtent( N, clifford_mat )
%RANDOMUNITARYSETWITHEXTENT generate set of random single-qubit unitaries,
%and compute their unitary extent.
% 
% Inputs:
%       N:  number of unitary operations to sample.
%
%       clifford_mat: matrix of choi state vectors for the single-qubit
%           Cliffords.
%
% unitary_set.U  standard matrix representations; complex 2-by-2-by-N array.
% unitary_set.PL  Pauli-Liouville representations; real 4-by-4-by-N array.
% unitary_set.SO  SO(3) representations; real 3-by-3-by-N array.
% unitary_set.extent  vector of unitary extent values;
% unitary_set.decomp  sum-over-Clifford decomposition, expressed as vector of coefficients; M-by-N array, where M is number of n-qubit Cliffords.
% unitary_set.status  solve status for the unitary extent problem.

[~,cliff_cols] = size(clifford_mat);


U_arr = zeros(2,2,N);
PL_arr = zeros(4,4,N);
SO_arr = zeros(3,3,N);
decomp_arr = zeros(cliff_cols,N);
extent_arr = zeros(N,1);
status_cell = cell(N,1);



for kk=1:N
    U_k = randU(2);
    U_arr(:,:,kk) = U_k;
    PL_k = PLdecomp(U_k,'unitary');
    PL_arr(:,:,kk) = PL_k;
    SO_k = PL_k(2:4,2:4);
    SO_arr(:,:,kk) = SO_k;
    [ extent_k,status_k,decomp_k,~,...
    ~] = findUnitaryExtent(U_k,clifford_mat,'high','SDPT3');
    decomp_arr(:,kk) = decomp_k;
    extent_arr(kk) = extent_k;
    status_cell(kk) = {status_k};
end

unitary_set.U = U_arr;
unitary_set.PL = PL_arr;
unitary_set.SO = SO_arr;
unitary_set.extent = extent_arr;
unitary_set.decomp = decomp_arr;
unitary_set.status = status_cell;


end

