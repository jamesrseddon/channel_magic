function [ state_vector ] = ketFromDensity( rho, varargin)
%KETFROMDENSITY Convert density matrix to state vector

if nargin > 1
    tolerance = varargin;
else
    tolerance =100*eps;
end

[rows,columns] = size(rho);

if rows ~= columns
    errorStruct.message = ['Input is not a square matrix.'...
                               ' Expecting a density matrix.'];
	errorStruct.identifier = ['quasi:ketFromDensity:'...
                                        'notSquare'];
	error(errorStruct);
end

rho_rank = rank(rho,tolerance);

if rho_rank ~= 1
    errorStruct.message = ['Input matrix is not rank 1.'...
                                'Expecting a pure state density matrix.'];
	errorStruct.identifier = ['quasi:ketFromDensity:'...
                                        'notPure'];
	error(errorStruct);    
end

[evecs,eval_mat] = eig(rho);

total_trace = trace(rho); % Checking the trace lets us find the state vector
                         % in the case when rho is not normalised.

[pure_index, ~] = find(abs(eval_mat - total_trace) <= tolerance);

state_vector = evecs(:,pure_index);

end

