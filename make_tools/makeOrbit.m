function [orbit_set] = makeOrbit(A,transform_set,varargin)
%MAKEORBIT make orbit of A under transformations in transform_set
%   Given a set of rotations transform_set and input matrix A, generate the
%   set whose elements are A_{j,k} = R_j*A*R_k, where R_i are the elements
%   of transform_set. varargin is an optional tolerance for checking
%   equivalence and removing duplicates. If second optional argument is
%   'P', elements equivalent up to global phase will also be removed.

if nargin > 2
    tolerance = varargin{1}
else
    tolerance = 100*eps
end

phase_flag = false;

if nargin > 3
    phase_option = varargin{2};
    phase_flag = strcmp(phase_option,'P')
end

num_transforms = size(transform_set,3);

next_index = 1;

orbit_set = [];

for jj = 1:num_transforms
    for kk = 1:num_transforms
        candidate = transform_set(:,:,jj)*A*transform_set(:,:,kk);
        if phase_flag
            duplicate_flag = checkMemberPhase(candidate,orbit_set);
        else
            duplicate_flag = checkMember(candidate,orbit_set,tolerance);
        end
        if ~duplicate_flag
            orbit_set(:,:,next_index) = candidate;
            next_index = next_index + 1;
        end
    end
end

total_unique = size(orbit_set,3);
disp(['Found ' num2str(total_unique) ' unique elements in the orbit.']);

end

