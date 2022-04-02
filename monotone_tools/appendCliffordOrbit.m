function [ new_unitary_set] = appendCliffordOrbit(  unitary_set, ...
    U, clifford_mat, clifford_cell,tol)
%APPENDCLIFFORDORBIT Append further gates to set of unitaries
% for upper bounding mixed channel extent. The appended gates are the 
% Clifford orbit of the input unitary U. Clifford_mat is the matrix of
% Clifford choi state vectors (needed for calculating extent of U.
% clifford_cell is cell array where column 2 gives Clifford gates as
% unitary matrices. tol is tolerance checking for duplicates in orbit.

% Calculate extent of U.

[ extent,~,~,~,~] = findUnitaryExtent(U,clifford_mat,'high','SDPT3');

clifford_arr = Ucell2array(clifford_cell);
tol
orbit_arr = makeOrbit(U,clifford_arr,tol,'P');

orbit_cell = Uarray2cell(orbit_arr);

new_unitary_set = appendUnitaryGates( unitary_set, ...
    orbit_cell, clifford_mat, extent);



end

