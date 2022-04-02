function [ pauli_vector, max_error ] = dyadPauliVector(ket,bra,pauli_array,varargin)
%DYADPAULIVECTOR Calculate Pauli vector for dyad |ket><bra|
%   Calculate the vector of expectation values:
%           <P_j> = Tr(P_j|ket><bra|) = <bra|P_j|ket>
%
% Values less than zero_tol will be zeroed. zero_tol is set to
% 10*eps by default, can be set manually by passing a number into
% varargin.

% check for varargin, check numeric and set zero_tol
if nargin > 3
    if isnumeric(varargin{1})
        zero_tol = varargin{1};
    else
        error_struct.message = ['Tolerance for killing off tiny values '...
            'should be numeric.'];
        error_struct.identifier = ['quasi:dyadPauliVector:'...
                'zeroToleranceNotNumeric'];
        error(error_struct);
    end
else
    zero_tol = 10*eps;
end

if ~(isvector(ket) && isvector(bra))
    error_struct.message = ['Input ket and bra states must be vectors.'];
    error_struct.identifier = ['channel_magic:dyadPauliVector:notVectors'];
    error(error_struct)
end

if length(ket) ~= length(bra)
    error_struct.message = ['Input ket and bra vectors should have same'...
                            ' dimension.'];
    error_struct.identifier = ['channel_magic:dyadPauliVector:'...
            'dimMismatch'];
    error(error_struct)
end

if ~(iscolumn(ket) && iscolumn(bra))
    error_struct.message = ['States to be used for bra and ket should '...
                            'be input as column vectors.'];
    error_struct.identifier = ['channel_magic:dyadPauliVector:'...
            'notColumns'];
    error(error_struct)
end

num_paulis = size(pauli_array,3);

raw_pauli_vector = zeros(num_paulis,1);

bra = bra';

for kk = 1:num_paulis
    this_pauli = pauli_array(:,:,kk);
    raw_pauli_vector(kk,1) = bra*this_pauli*ket;
end

real_vector = real(raw_pauli_vector);
imag_vector = imag(raw_pauli_vector);

clean_Im = imag_vector .* (abs(imag_vector) >= zero_tol);
clean_Re = real_vector .* (abs(real_vector) >= zero_tol);

pauli_vector = clean_Re + 1i*clean_Im;

max_error = max(abs(raw_pauli_vector - pauli_vector));

end

