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
function [ kraus ] = constructChannel(channel_type,params)
%CONSTRUCTCHANNEL Construct Kraus operators for a channel given a string to
% identify channel type and any parameters.
% Available configurations:
%   channel_type    params    description
%   'AMP_X'         p          Amplitude damping with noise parameter p, 
%                   theta       followed by X-rotation with angle theta
%   'X_AMP'         p          X-rotation followed by amplitude damping.
%                   theta
%   'Z_tensor'      n          n-fold tensor product of single-qubit
%                   theta       Z-rotations by angle theta.
% 
 
allowed_configs = {'AMP_X','X_AMP','Z_tensor'};

if ~any(strcmp(allowed_configs,channel_type))
    errorStruct.message = [channel_type ' is not an allowed channel type.'...
                ' Allowed strings are: ' sprintf('%s, ',allowed_configs{:})...
                ];
    errorStruct.identifier = ['quasi:constructChannel:invalidChannelType'];
    error(errorStruct)
end

switch channel_type
    case 'AMP_X'
        if ~(isfield(params,'p') && isfield(params,'theta'))
            errorStruct.message = ['The channel specified requires '...
                'parameters to include noise ''p'' and angle ''theta''.'];
            errorStruct.identifier = ['quasi:constructChannel:'...
                'missingParameters'];
            error(errorStruct);
        end
        p = params.p;
        theta = params.theta;
        display(['Constructing amplitude damping with noise p = '...
            num2str(p)...
            ', followed by X-rotation through angle ' num2str(theta)'.']);
        
        amp_damp = amplitudeDamping(p);
        X_rotate = singleQubitPauliRot(theta,'X');
        kraus = composeChannels(amp_damp,X_rotate);
    case 'X_AMP'
        if ~(isfield(params,'p') && isfield(params,'theta'))
            errorStruct.message = ['The channel specified requires '...
                'parameters to include noise ''p'' and angle ''theta''.'];
            errorStruct.identifier = ['quasi:constructChannel:'...
                'missingParameters'];
            error(errorStruct);
        end
        p = params.p;
        theta = params.theta;
        display(['Constructing X-rotation through angle ' num2str(theta)...
            ', followed by amplitude damping with noise p = ' num2str(p)...
             '.']);
        
        amp_damp = amplitudeDamping(p);
        X_rotate = singleQubitPauliRot(theta,'X');
        kraus = composeChannels(X_rotate,amp_damp);
    case 'Z_tensor'
        if ~(isfield(params,'n') && isfield(params,'theta'))
            errorStruct.message = ['The channel specified requires '...
                'parameters to include number of qubits ''n'' and angle'...
                ' ''theta''.'];
            errorStruct.identifier = ['quasi:constructChannel:'...
                'missingParameters'];
            error(errorStruct);
        end
        n = params.n;
        theta = params.theta;
        display(['Constructing Z-rotation through angle ' num2str(theta)...
            'with number of qubits n = '...
            num2str(n)...
            '.']);
        
        Z_rotate = singleQubitPauliRot(theta,'Z');
        current_Z = Z_rotate;
        for kk = 1:(n-1)
            current_Z = kron(current_Z,Z_rotate);
        end
        kraus = current_Z;
end