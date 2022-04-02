function [ SO_mat, E_vec, status, Xi_bound, witness, optbnd, delta_vec, xi_vec]...
    = upperBoundChannelExtent(channel, unitary_set)
%UPPERBOUNDCHANNELEXTENT Upper bound mixed channel extent for single-qubit
%channels, using finite subset of witness constraints.
%   The dual formulation for mixed channel extent Xi(M) can be written:
%       minimise    Tr[W M]
%       subject to  xi(R) >= Tr[W R] for all R in SO(3)
%   where M is the Pauli transfer matrix for the single-qubit channel of
%   interest, R are Pauli transfer matrices for SO(3) rotations in the
%   Pauli basis, and xi is the unitary extent. We can't solve directly as
%   SO(3) is an infinite set, but we can upper bound Xi(M) by reducing to a
%   finite number of constraints based on a finite subset of SO(3).
%   Inputs:
%     channel:      target single-qubit channel (currently assumed to be
%                   depolarised single-qubit rotation, in Pauli-Liouville
%                   representation).
%     unitary_set:  Subset of SO(3) rotations, with pre-computed unitary
%                   extent. unitary_set should be a structure with fields:
%                       unitary_set.SO (SO(3) representation; 3-by-3-by-M 
%                       array, where M is number of rotations);
%                       unitary_set.extent (xi(R) calculated for each
%                                           rotation R).
%
delta = 0

if isfield(unitary_set,'SO')
    SO_arr = unitary_set.SO;
else
    error_struct.message = ['SO field missing from input structure.'];
    error_struct.identifier = 'quasi:upperBoundChannelExtent:SOmissing';
    error(error_struct);
end

if isfield(unitary_set,'extent')
    xi_vec = unitary_set.extent;
else
    error_struct.message = ['Extent field missing from input structure.'];
    error_struct.identifier = 'quasi:upperBoundChannelExtent:extentMissing';
    error(error_struct);
end

[channel_rows, channel_cols] = size(channel);
[SO_rows, SO_cols, SO_num] = size(SO_arr);
[extent_num, extent_cols] = size(xi_vec);

if ~(channel_rows == 4 && channel_cols == 4)
    error_struct.message = ['Input channel should be represented as '...
                                '4-by-4 Pauli transfer matrix.'];
    error_struct.identifier = 'quasi:upperBoundChannelExtent:invalidChannel';
    error(error_struct);
end

if ~(SO_rows == 3 && SO_cols == 3)
    error_struct.message = ['SO field should be 3-by-3-by-M array.'];
    error_struct.identifier = 'quasi:upperBoundChannelExtent:invalidSO';
    error(error_struct);
end

if extent_cols ~= 1
    error_struct.message = ['Extent field should be a vector.'];
    error_struct.identifier = 'quasi:upperBoundChannelExtent:invalidExtent';
    error(error_struct);
end

if extent_num ~= SO_num
    error_struct.message = ['Number of SO(3) matrices and number of '...
                            'extent values do not match.'];
    error_struct.identifier = 'quasi:upperBoundChannelExtent:countMismatch';
    error(error_struct);
end

% Reshape the SO(3) array so that it is a matrix where each column is a
% 10-element vector. The 1 in the first row ensures correct channel
% normalisation in the PL representation.

SO_mat = zeros(10,SO_num);

SO_mat(1,:) = 1; 

for kk = 1:SO_num
    current_SO = SO_arr(:,:,kk);
    SO_mat(2:10,kk) = reshape(current_SO,9,1);
end

SO_mat = SO_mat';

% Similarly reshape the PL transfer matrix for the target channel.

E_vec = reshape(channel,16,1);
E_vec = E_vec([1 6:8 10:12 14:16]); % Discard zero entries in first row and 
                    % column of PL matrix, as we assumed that it represents
                    % a depolarised single-qubit rotation.
E_vec = E_vec';
     
delta_vec = ones(SO_num,1)*delta;
                    
cvx_begin 
    variable q(10);
    maximize( E_vec*q );
    subject to
        SO_mat * q <= xi_vec + delta_vec
cvx_end

status = cvx_status;

if strcmp(status,'Solved')||strcmp(status,'Inaccurate/Solved')
    Xi_bound = cvx_optval;
    witness = q;
    dual_var = 0;
    if exist('cvx_optbnd');
        optbnd = cvx_optbnd;
    elseif exist('cvx_bound')
        optbnd = cvx_bound;
    end
else
    Xi_bound = 0;
    witness = q;
    dual_var = 0;
    optbnd = 0;
end

