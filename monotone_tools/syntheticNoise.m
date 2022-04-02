function [dyadic_negativity,rhoprime,status,dual_bound] = ...
    syntheticNoise(rho,omega,lambda,A_matrix)
% SYNTHETICNOISE Calculate a noisy approximation to a state then find dyadic negativity.
%	Given a state of interest rho, some other state omega, and a mixing 
%   parameter lambda, calculate the mixture rho' of rho with omega. Then
%   calculate the dyadic negativity of rho'.

rhoprime = (rho + (lambda - 1)*omega)/lambda;

[dyadic_negativity,status,distribution,dual_var,...
    dual_bound] = findRobustness(rhoprime,A_matrix,'high','SDPT3','C');

