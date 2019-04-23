% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
%  Copyright © 2019 James R. Seddon
%  This file is part of project: Channel Magic
%
function [ L_of_rho ] = unitaryGenerator( input_state,hamiltonian )
%UNITARYGENERATOR Given a Hamiltonian, apply the corresponding generator of
% unitary evolution to an input_state
%   Quantum states evolve according to the Liouville-Von Neumann equation:
%    d(RHO)/dt = -i [ H, RHO] = L(RHO), where H is the system Hamiltonian,
%   and L is the superoperator defined by:
%   L(A) = -i [ H, A] 
%   For time-independent L, the evolution equation has formal solution:
%   RHO(t) = exp(L t) [RHO(0)]
%   We call L the generator of the evolution. (This can be generalised for
%   non-unitary processes.)
%   This function takes as input H and RHO and calculates L(RHO).

% Check inputs Hermitian up to some tolerance.
herm_tol = 10*eps;

if ~checkNearlyHermitian(input_state,herm_tol)
    error_struct.message = ['Input state was not Hermitian.'];
    error_struct.identifier = 'quasi:unitaryGenerator:stateNotHermitian';
    error(error_struct);
end

if ~checkNearlyHermitian(hamiltonian,herm_tol)
    error_struct.message = ['Hamiltonian was not Hermitian.'];
    error_struct.identifier = ['quasi:unitaryGenerator:'...
                                    'hamiltonianNotHermitian'];
    error(error_struct);
end

L_of_rho = -1i*(hamiltonian*input_state - input_state*hamiltonian);


end

