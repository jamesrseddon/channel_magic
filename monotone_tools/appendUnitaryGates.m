function [ new_unitary_set ] = appendUnitaryGates( unitary_set, ...
    unitary_cell, clifford_mat, varargin)
%APPENDUNITARYGATES Append further unitary gates to set of unitaries
% unitaries for upper bounding mixed channel extent. Optionally pass in
% vector of extent values (or single value) if these are already known, 
% to avoid having to recompute. N.B. in this case the decompositions will
% not be computed.

known_extent = false;

extent_input = 0;

if nargin > 3
    known_extent = true;
    extent_input = varargin{1};
end

U_orig = unitary_set.U;

M = size(U_orig,3);

M_app = size(unitary_cell,1);

if known_extent && length(extent_input) == 1
    extent_input = extent_input*ones(M_app,1);
end

M_Cl = size(clifford_mat,2);

N = M + M_app;

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



for kk=1:M_app
    new_kk = M + kk;
    U_k = unitary_cell{kk,2};
    U_arr(:,:,new_kk) = U_k;
    PL_k = PLdecomp(U_k,'unitary');
    PL_arr(:,:,new_kk) = PL_k;
    SO_k = PL_k(2:4,2:4);
    SO_arr(:,:,new_kk) = SO_k;
    if known_extent
        extent_arr(new_kk) = extent_input(kk);
        status_cell{new_kk} = 'precomputed';
    else
       [ extent_k,status_k,decomp_k,~,...
                ~] = findUnitaryExtent(U_k,clifford_mat,'high','SDPT3');
        extent_arr(new_kk) = extent_k;    
        decomp_arr(:,new_kk) = decomp_k;
        status_cell(new_kk) = {status_k};
    end
end

new_unitary_set.U = U_arr;
new_unitary_set.PL = PL_arr;
new_unitary_set.SO = SO_arr;
new_unitary_set.extent = extent_arr;
new_unitary_set.decomp = decomp_arr;
new_unitary_set.status = status_cell;

end

