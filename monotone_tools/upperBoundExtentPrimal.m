function [extent_bound,status,distribution,dual_var,...
    optbnd] = upperBoundExtentPrimal(state,A_matrix,extent,varargin)
%UPPERBOUNDEXTENTPRIMAL Compute upper bound to mixed extent using primal
%formulation.
%   Optimisation problem:
%       variable p
%       minimise  p' * E
%       subject to A*diag(p)*A'
%                   p > = 0
%       
%       where A is a matrix s. t. each column is state vector with known
%       extent. E is the vector of these extent values.
%   'state' should be passed in as density matrix in standard basis.
%   'A_matrix' is list of states with known extent.
%   If an argument 'high' is passed in, CVX will aim for higher precision.
%   Second optional argument specifies the solver to be used by CVX. Pass
%   in string 'SDPT3' to use that solver, otherwise SeDuMi will be used.

[A_rows, N] = size(A_matrix);
[state_rows, state_cols] = size(state)
num_extent = length(extent);

high_precision = false;
flags = varargin(:);
num_flags = length(flags);
% Sets a flag for high_precision in solver if 'high' is the third argument.
% class(flags)
if num_flags > 0
    if strcmp(flags{1},'high')
        high_precision = true;
    end
end

SDPT_flag = false;
if num_flags > 1
    if strcmp(flags{2},'SDPT3')
        SDPT_flag = true;
    end
end

if state_rows ~= state_cols
    error_struct.message = ['First argument should be a square matrix.'];
    error_struct.identifier = ['quasi:upperBoundExtentPrimal:' ...
                                    'unSquareDensityMatrix'];
    error(error_struct);
end

if state_rows ~= A_rows
    error_struct.message = ['State vector matrix dimension does not '...
                'match target state.'];
    error_struct.identifier = ['quasi:upperBoundExtentPrimal:' ...
                                    'dimMismatch'];
    error(error_struct);
end


if num_extent ~= N
    error_struct.message = ['Expected ' num2str(N) ' but extent vector '...
        'has ' num2str(num_extent) ' entries.'];
    error_struct.identifier = ['quasi:upperBoundExtentPrimal:'...
        'numberMismatch'];
    error(error_struct);
end

    cvx_begin
            cvx_solver sedumi
            variable p(N)
            P = diag(p)
            dual variable y
            minimize(p'*extent)
            subject to
                y: A_matrix*P*A_matrix' == state
                   p >= 0
    cvx_end

 %%% run CVX linear program.

solve_type = 'sedumi_low';

if high_precision && ~SDPT_flag
    solve_type = 'sedumi_high';
end

if ~high_precision && SDPT_flag
    solve_type = 'SDPT_low';
end

if high_precision && SDPT_flag
    solve_type = 'SDPT_high';
end


display(['Solve type is ' solve_type '.']);

switch solve_type
    case 'sedumi_low'
        cvx_begin
            cvx_solver sedumi
            variable p(N)
            P = diag(p)
            dual variable y
            minimize(p'*extent)
            subject to
                y: A_matrix*P*A_matrix' == state
                   p >= 0
        cvx_end
        
    case 'sedumi_high'
        
        cvx_begin
            cvx_solver sedumi
            cvx_precision high
            variable p(N)
            P = diag(p)
            dual variable y
            minimize(p'*extent)
            subject to
                y: A_matrix*P*A_matrix' == state
                   p >= 0
        cvx_end
    case 'SDPT_low'
        cvx_begin
            variable p(N)
            P = diag(p)
            dual variable y
            minimize(p'*extent)
            subject to
                y: A_matrix*P*A_matrix' == state
                   p >= 0
        cvx_end
    case 'SDPT_high'
        cvx_begin
            cvx_precision high
            variable p(N)
            P = diag(p)
            dual variable y
            minimize(p'*extent)
            subject to
                y: A_matrix*P*A_matrix' == state
                   p >= 0
        cvx_end
end


status = cvx_status;

if strcmp(status,'Solved')||strcmp(status,'Inaccurate/Solved')
    extent_bound = cvx_optval;
    distribution = p;
    dual_var = y;
    if exist('cvx_optbnd');
        optbnd = cvx_optbnd;
    elseif exist('cvx_bound')
        optbnd = cvx_bound;
    end
else
    extent_bound = 0;
    distribution = 0;
    dual_var = 0;
    optbnd = 0;
end 
    
    
end

