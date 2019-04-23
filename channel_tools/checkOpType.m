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
function op_type = checkOpType( input_array )
%CHECKOPTYPE Decide whether array represents unitary or Kraus operators.
%   Takes as input an array, which can be 2D or 3D. If 2D, checks that it
%   is a valid unitary. If 3D, checks whether forms a complete Kraus
%   representation (third dimension labels the Kraus operator).

tolerance = 50*eps; % tolerance for checking unitarity/completeness.
num_dims = ndims(input_array);

if (num_dims ~= 2) && (num_dims ~= 3)
    errorStruct.message = [ 'Input array has ' num2str(num_dims)...
        ' dimensions, but checkOpType cannot interpret arrays with '...
        'more than three dimensions.'];
    errorStruct.identifier = ['quasi:checkOpType:'...
        'unrecognisedDimension'];
    error(errorStruct);
end

[rows,columns,~] = size(input_array);

if rows ~= columns
    errorStruct.message = ['Expect either a square matrix or a 3D array'...
        ' where first two dimensions are equal, but input has '...
        num2str(rows) ' rows and ' num2str(columns) ' columns.'];
    errorStruct.identifier = ['quasi:checkOpType:'...
        'notSquare'];
    error(errorStruct);
end

identity = eye(rows);

switch num_dims
    case 2
        % Check for unitarity.
        UdaggerU = input_array'*input_array;
        is_unitary = checkEqual(UdaggerU,identity,tolerance);
        if is_unitary
            display('Input is a unitary matrix.');
            op_type = 'U';
        else
            display('Input is non-unitary matrix.');
            op_type = 'NU';
        end
    case 3
        % Check for complete Kraus representation.
        num_operators = size(input_array,3);
        for kk = 1:num_operators
            K = input_array(:,:,kk);
            kraus_product(:,:,kk) = K'*K;
        end
        kraus_sum = sum(kraus_product,3);
        is_complete = checkEqual(kraus_sum,identity,tolerance);
        if is_complete
            display('Input forms a complete Kraus representation.');
            op_type = 'K';
        else
            display('Input does not form complete Kraus representation.');
            op_type = 'NK';
        end
end

end
