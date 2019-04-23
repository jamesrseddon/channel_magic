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
function [output] = pauliTensor(paulistring)
%%%% Forms the matrix corresponding to the tensor product of n Pauli
%%%% operators based on an input string using standard Pauli notation,
%%%% Example:
%%%% string = 'Z X Z I I I'
%%%% corresponds to Z on qubit 1, X on 2, Z on 3 and I on remaining.


parsedString = textscan(paulistring,'%s');

pauliList = parsedString{1};

qubitnumber = length(pauliList);

% Define scalars
minus = -1;
plusi = 1i;
minusi = -1i;

% Define the Paulis

X = [ 0 1;
      1 0];
Y = 1i*[0 -1;
        1 0];
Z = [1 0;
     0 -1];
I = eye(2);

% Define projectors onto Pauli eigenstates
X0 = (1/2)*[1 1; 
                 1 1];
X1 = (1/2)*[1 -1;
                  -1 1];
              
Y0 = (1/2)*[1 -1i;
                  1i 1];
Y1 = (1/2)*[1 1i;
                  -1i 1];

Z0 = [1 0;
     0 0];
Z1 = [0 0;
      0 1];


% Create a key-value map

pauliMap = containers.Map();
pauliMap('I') = I;
pauliMap('X') = X;
pauliMap('Y') = Y;
pauliMap('Z') = Z;

pauliMap('X0') = X0;
pauliMap('X1') = X1;
pauliMap('Y0') = Y0;
pauliMap('Y1') = Y1;
pauliMap('Z0') = Z0;
pauliMap('Z1') = Z1;

pauliMap('-') = minus;
pauliMap('i') = plusi;
pauliMap('-i') = minusi;

firstPauli = pauliMap(pauliList{1});

currentTensor = firstPauli;

for k = 2:qubitnumber
    currentTensor = kron(currentTensor,pauliMap(pauliList{k}));
end

output = currentTensor;