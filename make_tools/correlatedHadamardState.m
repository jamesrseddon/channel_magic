function [ RHO ] = correlatedHadamardState( p )
%CORRELATEDHADAMARDSTATE constructs state: 
% RHO = p |H0 H0><H0 H0| + (1-p)|H1 H1><H1 H1|
%
%   H0 and H1 are the eigenvectors of the Hadamard.

if p < 0 || p > 1 
    error_struct.message = ['Argument must be in the range 0 <= p <= 1'];
    error_struct.identifier = ['channel_magic:correlatedHadamardState:'...
                                    'invalidProbability'];
    error(error_struct)
end

Hadamard = 1/sqrt(2)*[1 1; 1 -1];

[V, D] = eigs(Hadamard);

H0 = V(:,2);

H1 = V(:,1);


H0H0 = kron(H0,H0);
H1H1 = kron(H1,H1);

RHO = p*H0H0*H0H0' + (1-p)*H1H1*H1H1';

end

