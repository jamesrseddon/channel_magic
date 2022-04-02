function [lambda,x,sigma,omega] = RobMagGdecomp(rho,stabiliser_matrix)
%ROBMAGGDECOMP computes the generalised robustness of magic and gives
%associated state decomposition.
%   requires as second argument stabiliser_matrix whose columns are pure 
%   stabilizer states. Computes generalised robustness lambda, then gives
%   the density matrices making up the associated decomposition:
%       omega = (lambda*sigma - rho)/(lambda - 1)
%
%   Requires function RobMagG rom Bartosz Regula's code accompanying 
%   PRX Quantum 2, 010345 (2021), available at:
%   https://bartoszregula.me/code/magic (accessed 2nd April 2022)

[lambda, x] = RobMagG(rho,stabiliser_matrix);

sigma = stabiliser_matrix*diag(x)*stabiliser_matrix'/lambda;

omega = (lambda*sigma - rho)/(lambda - 1);

end

