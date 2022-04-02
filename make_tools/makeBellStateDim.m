function [ bell_state] = makeBellStateDim(d)
%MAKEBELLSTATEDIM Construct Bell state of local dimension d.
%   Make the state vector in standard basis for maximally entangled state:
%       |Phi> = (1/sqrt(d)) \sum_{j=1}^d |j> |j>
%   i.e. the local dimension is d, dimension of the entangled state is d^2.

if ~(mod(d,1) == 0 && d > 0 && length(d) == 1)
    error_struct.message = ['Input dimension not a natural number'];
    error_struct.identifier = ['quasi:makeBellStateDim:invalidDimension'];
    error(error_struct)
end

working_vector = zeros(d,1);

for kk = 1:d
    current_index = (kk-1)*d + kk;
    working_vector(current_index) = 1;
end

bell_state = (1/sqrt(d))*working_vector;

end


