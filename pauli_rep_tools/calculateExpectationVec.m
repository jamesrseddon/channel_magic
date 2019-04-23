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
function [output,real_vector] = calculateExpectationVec(state,observables,varargin)
%%% Takes as input a state and an array of observables, and outputs a
%%% vector of expectation values. State should be in density matrix
%%% format. Observables should be a 3-dimensional array, where the third
%%% dimension is the index of the observable, eg:
%%%
%%% observables(:,:,k) = O_k
%%%
%%% where O_k is a matrix with same dimension as the input state.
%%% As a result of numerical error, calculated expectation values can have
%%% small imaginary components - these will be cleaned up. It is also
%%% possible to get tiny non-zero values when the true value should be
%%% zero. This can make optimisation problems appear to be infeasible when
%%% in fact they should not be.
%%% Values less than zero_tol will be zeroed. zero_tol is set to
%%% 10*eps by default, can be set manually by passing a number into
%%% varargin.

% check for varargin, check numeric and set zero_tol
if nargin > 2
    if isnumeric(varargin{1})
        zero_tol = varargin{1};
    else
        error_struct.message = ['Tolerance for killing off tiny values '...
            'should be numeric.'];
        error_struct.identifier = ['quasi:calculateExpectationVec:'...
                'zeroToleranceNotNumeric'];
        error(error_struct);
    end
else
    zero_tol = 10*eps;
end

% check input variable state is a valid density matrix
if ~(ndims(state)== 2)
    errorStruct.message = ['State needs to be input as density matrix '...
                    'but input array was not a matrix.'];
    errorStruct.identifier = ['quasi:calculateExpectationVec:'...
                'inputStateNotMatrix'];
    error(errorStruct);
end

[rows columns] = size(state);

% check square
if rows ~= columns
    errorStruct.message = ['State needs to be input as density matrix '...
                    'but input matrix was not square.'];
    errorStruct.identifier = ['quasi:calculateExpectationVec:'...
                'inputStateNotSquare'];
    error(errorStruct);
end

% check Hermitian
if ~checkEqual(state,state',1e-14)
    errorStruct.message = ['State needs to be input as density matrix '...
                    'but input matrix was not hermitian.'];
    errorStruct.identifier = ['quasi:calculateExpectationVec:'...
                'inputStateNotHermitian'];
    error(errorStruct);
end

% check observables array has right dimensionality.
if ~(ndims(observables)== 3)
    errorStruct.message = ['Observables need to be input as '...
                    'three-dimensional array, where third index labels '...
                    'each observable, but input array was not 3D.'];
    errorStruct.identifier = ['quasi:calculateExpectationVec:'...
                'inputObservablesNot3D'];
    error(errorStruct);
end

[obs_rows, obs_columns, num_observables] = size(observables);

% check observables array has right dimensionality.
if ~((obs_rows == rows) && (obs_columns == columns))
    errorStruct.message = ['Observables did not have same dimension as '...
                                'the input state.'];
    errorStruct.identifier = ['quasi:calculateExpectationVec:'...
                'inputObservablesDimensionMismatch'];
    error(errorStruct);
end

exp_vector = zeros(num_observables,1);


for oo = 1:num_observables;
    current_obs = observables(:,:,oo);
    exp_vector(oo,1) = trace(current_obs*state);
end

real_vector = real(exp_vector);

clean_vector = real_vector .* (abs(real_vector) >= zero_tol);

clean_vector - real_vector;

output = clean_vector;

end

    