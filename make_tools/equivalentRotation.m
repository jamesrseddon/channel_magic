function [ phi , rotation_matrix, U, extent] = equivalentRotation( alpha,p)
%EQUIVALENTROTATION 

phi = (asin((1-p)*sin(2*alpha + pi/4) + p/(2*sqrt(2))) - pi/4)/2;


rotation_matrix = [ cos(2*phi) -sin(2*phi) 0;
                    sin(2*phi) cos(2*phi) 0;
                    0 0 1];

 U = [exp(-1i*phi) 0 ;
      0 exp(1i*phi)];

extent = (cos(pi/4 - phi) + (sqrt(2)-1)*sin(pi/4 - phi))^2;
end

