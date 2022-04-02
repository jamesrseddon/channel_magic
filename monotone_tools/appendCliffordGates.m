function [ new_unitary_set ] = appendCliffordGates( unitary_set, unitary_cell)
%APPENDCLIFFORDGATES Append Clifford gates to set of unitaries
% for upper bounding mixed channel extent. unitary_cell is the
% cell array containing Cliffords in column 2.

U_orig = unitary_set.U;

M = size(U_orig,3);

M_Cl = size(unitary_cell,1);

N = M + M_Cl;

U_arr = zeros(2,2,N);
PL_arr = zeros(4,4,N);
SO_arr = zeros(3,3,N);
decomp_arr = zeros(M_Cl,N);
extent_arr = zeros(N,1);
status_cell = cell(N,1);

U_arr(:,:,1:M) = U_orig;
PL_arr(:,:,1:M) = unitary_set.PL;
SO_arr(:,:,1:M) = unitary_set.SO;
decomp_arr(:,1:M) = unitary_set.decomp;
extent_arr(1:M) = unitary_set.extent;
status_cell(1:M) = unitary_set.status;



for kk=1:M_Cl
    new_kk = M + kk;
    U_k = unitary_cell{kk,2};
    U_arr(:,:,new_kk) = U_k;
    PL_k = PLdecomp(U_k,'unitary');
    PL_arr(:,:,new_kk) = PL_k;
    SO_k = PL_k(2:4,2:4);
    SO_arr(:,:,new_kk) = SO_k;
    decomp_arr(kk,new_kk) = 1;
    extent_arr(new_kk) = 1;
    status_cell(new_kk) = {'clifford'};
end

new_unitary_set.U = U_arr;
new_unitary_set.PL = PL_arr;
new_unitary_set.SO = SO_arr;
new_unitary_set.extent = extent_arr;
new_unitary_set.decomp = decomp_arr;
new_unitary_set.status = status_cell;

end

