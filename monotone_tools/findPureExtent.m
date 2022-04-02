function [ extent,status,distribution,dual_var,...
    optbnd] = findPureExtent( state_vector, A_matrix, varargin )
%FINDPUREEXTENT Find the stabiliser extent for a given pure state vector.
%   Given a state |psi>, find min{(||c||_1)^2 : |psi> = \sum_j c_j |phi_j>}
%   where |phi_j> are stabiliser states. |psi> should be represented as a
%   vector in the computational basis. A_matrix should be a matrix where
%   each column is a state vector for one of the pure stabiliser states, of
%   the same dimension as the state of interest.
%   If an argument 'high' is passed in, CVX will aim for higher precision.
%   Second optional argument specifies the solver to be used by CVX. Pass
%   in string 'SDPT3' to use that solver, otherwise SeDuMi will be used.
%   For pure state extent, complex variables are always enabled.

dim_input = size(state_vector,1);

num_qubits = log2(dim_input);
[A_rows A_columns] = size(A_matrix);

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

complex_flag = true;

if dim_input ~= A_rows
    errorStruct.message = ['Expected ' num2str(dim_input) ' but specified A matrix '...
        'has ' num2str(A_rows) ' rows.'];
    errorStruct.identifier = 'quasi:findPureExtent:mismatchedAmatrix';
    error(errorStruct);
end

%%% calculate b vector for the input state.

b_vector = state_vector;


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

if complex_flag
    solve_type = [solve_type '_complex'];
end

display(['Solve type is ' solve_type '.']);

switch solve_type
    case 'sedumi_low_complex'
        cvx_begin
            cvx_solver sedumi
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
        
    case 'sedumi_high_complex'
        
        cvx_begin
            cvx_solver sedumi
            cvx_precision high
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'SDPT_low_complex'
        cvx_begin
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
    case 'SDPT_high_complex'
        cvx_begin
            cvx_precision high
            variable q(A_columns) complex;
            dual variable y;
            minimize( norm(q,1) );
            subject to
                y : A_matrix * q == b_vector
        cvx_end
end


status = cvx_status;

if strcmp(status,'Solved')||strcmp(status,'Inaccurate/Solved')
    extent = cvx_optval^2;
    distribution = q;
    dual_var = y;
    if exist('cvx_optbnd');
        optbnd = cvx_optbnd^2;
    elseif exist('cvx_bound')
        optbnd = cvx_bound^2;
    end
else
    extent = 0;
    distribution = 0;
    dual_var = 0;
    optbnd = 0;
end


end

