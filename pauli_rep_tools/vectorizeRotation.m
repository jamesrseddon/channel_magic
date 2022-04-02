function [R_vec] = vectorizeRotation(R)
%VECTORIZEROTATION Vectorize an SO(3) rotation matrix.
%   

R_vec = reshape(R,9,1);

end

