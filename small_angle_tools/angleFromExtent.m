function [a, U_a, psi_a, R_a ] = angleFromExtent(extent)
%ANGLEFROMEXTENT Find a single-qubit Z-rotation with given extent.
%   Input is the pure state extent for the Choi state:
%   |a> = U_a |+>, where:
%   U_a = exp(-i*Z*a) = [ exp(-i*a) 0;
%                         0        exp(i*a)];
%   We then calculate the angle a, using the fact that:
%   extent(|a>) = (cos(pi/4 - a) + (sqrt(2) - 1)*sin(pi/4 - a))^2.
%   Output angle a, unitary U_a, state |a>, SO(3) representation for U_a.

c = 3*pi/8;

b = sqrt(2 + sqrt(2))/2;

a = asin(sqrt(extent)*b) - c;

U_a = [exp(-1i*a) 0;
       0 exp(1i*a)];

plus = (1/sqrt(2))*[1;1];

psi_a = U_a * plus;

R_a = [ cos(2*a) -sin(2*a) 0;
        sin(2*a) cos(2*a) 0;
        0 0 1];

end

