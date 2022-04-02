function [ A_matrix, label_cell, error_vector] = ...
                makeAmatrixDyads( dyad_set, num_qubits, varargin)
%MAKEAMATRIXDYADS Calculate A matrix for a set of stabiliser state dyads.
%   Input should be of the form:
%       {ket_vec_1, bra_vec_1, dyad_1_1, ket_stabiliser_1, bra_stabiliser_1;
%       ket_vec_1, bra_vec_2, dyad_1_2, ket_stabiliser_1, bra_stabiliser_2;
%            ...}
%
%   jth row of label_cell gives the stabiliser group for the ket (column 1) 
%   and bra (column 2) of the jth column of A_matrix.
%
%   error_vector lists the max element-wise absolute error found in
%   cleaning up tiny numerical errors.
%
%   varargin can be a tolerance for killing off tiny numerical errors.

if nargin > 2
    if isnumeric(varargin{1})
        zero_tol = varargin{1};
    else
        error_struct.message = ['Tolerance for killing off tiny values '...
            'should be numeric.'];
        error_struct.identifier = ['quasi:makeAmatrixDyads:'...
                'zeroToleranceNotNumeric'];
        error(error_struct);
    end
else
    zero_tol = 10*eps;
end

[pauli_array,~] = enumeratePaulis(num_qubits);

num_dyads = size(dyad_set,1);

num_paulis = 4^num_qubits;

A_matrix = zeros(num_paulis,num_dyads);
label_cell = dyad_set(:,4:5);
error_vector = zeros(num_dyads,1);

perc_complete = 0;
fprintf(1,'Building %d A matrix columns\n', num_dyads);
fprintf(1,'Columns calculated: %3d%%\n', perc_complete);

for kk = 1:num_dyads
    ket = dyad_set{kk,1};
    bra = dyad_set{kk,2};
    [this_vector, max_err] = dyadPauliVector(ket,bra,pauli_array,zero_tol);
    A_matrix(:,kk) = this_vector;
    error_vector(kk,1) = max_err;
    
    perc_complete = 100*kk/num_dyads;
    fprintf(1,'\b\b\b\b%3.0f%%', perc_complete);
end
fprintf(1,'\n');
end

