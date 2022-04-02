function [xi,alpha,T,Lambda,Lambda_plus,xi_theta,Q,Pr,status,A_orbit] = ...
    equimagicalSearch(channel_type,params,search_interval,tolerance,...
    A_comp_basis_1,A_comp_basis_2,A_dyad)
%EQUIMAGICALSEARCH Search for upper bound of mixed channel extent.
%   Given channel E, Search for an upper bound of Xi(E) by looking for
%   equimagical decompositions comprised of elements from the Clifford
%   orbit of U(alpha) = exp[-i*Z*alpha], where alpha is swept through a
%   range of values corresponding to intervals in the extent.
%   Inputs are:
%       channel_type = string to specify channel as per constructChannel,
%       eg. 'Z_DEP' for depolarised Z-rotation.
%       params = channel parameters as per constructChannel, eg. for
%       'Z_DEP', this is specified as:
%               params.theta = rotation angle in radians
%               params.p = noise parameter.
%
%       search_interval = interval for incrementing the pure state extent, 
%       xi(U(alpha)).
%       tolerance = termination criteria - algorithm terminates when
%       L1-norm of decomposition < 1 + tolerance.
%       
%       A_comp_basis_1 = Matrix where each column is the state vector for a
%       single-qubit stabiliser state.
%       A_comp_basis_2 = As per A_comp_basis_1, but for two-qubit states.
%       A_dyad = Matrix where each column is Pauli representation for
%       two-qubit stabiliser state dyad.
%
%       Initial extent xi_0 is the dyadic negativity E(Lambda). Then at
%       the t-th step the extent will be:
%           xi_t = xi_0 + search_interval*t.
%       The t-th alpha value is computed so that xi(U(alpha_t)) = xi_t.
%       The search terminates when an equimagical decomposition is found,
%       or when xi_t = xi(U(theta)) (i.e. the extent of the noiseless
%       rotation).
%       
%       Outputs:
%           xi = final extent (upper bounds Xi(E)).
%           alpha = final rotation angle.
%           T = number of iterations.
%           Lambda = dyadic negativity of E.
%           Lambda_plus = generalised robustness of E.
%           xi_theta = extent of noiseless rotation U(theta).
%           Q = l1-norm of decomposition in final step. If the search is 
%               successful this will be 1 (up to the specified tolerance).
%           Pr = Final probability distribution over elements in the
%               Clifford orbit (or quasiprobability if unsuccessful).
%           A_orbit = Matrix where each column is Pauli representation of 
%               elements in the orbit. 
%           status = CVX solve status of the final run of linear program.

[ kraus ] = constructChannel(channel_type,params);

y = params.theta;

U_y = [exp(-1i*y) 0;
       0 exp(1i*y)];

plus = (1/sqrt(2))*[1;1];

psi_y = U_y * plus;

[xi_theta,~,~,~,~] = findPureExtent( psi_y,...
                             A_comp_basis_1,'high','SDPT3','C');
choi_state = makeChoiState(kraus);
[Lambda_plus,~] = RobMagG(choi_state,A_comp_basis_2);
                                          


[Lambda,~,~,~,~] = ...
    findRobustness(choi_state,A_dyad,'high','SDPT3','C');

xi_0 = Lambda;
xi_vec = [];
alpha_vec = [];
Q_vec = [];
delta = search_interval;
t = 0;
prob_dist_tol = tolerance;

continue_search = true;

while continue_search
    xi_t = xi_0 + delta*t;
    xi_vec = [xi_vec;xi_t];
    [alpha_t, ~, ~,~ ] = angleFromExtent(xi_t);
    alpha_t = real(alpha_t);
    alpha_vec = [alpha_vec; alpha_t];
    Amat_t = makeAmatrixOrbit(alpha_t);
    [Q_t,status_t,dist_t,~,dist_optbnd] = ...
        findRobustness(choi_state,Amat_t,'high','SDPT3');
    Q_vec = [Q_vec;Q_t];
    if and(Q_t <= 1 + prob_dist_tol,Q_t > 0)
        continue_search = false;
        success = true;
    end
    if alpha_t > xi_theta
        continue_search = false;
        success = false;
    end
    t = t + 1;
end

T = length(Q_vec);

xi = xi_t;
alpha = alpha_t;
Q = Q_t;
Pr = dist_t;
A_orbit = Amat_t;
status = status_t;

if success
    disp(['Iterations: ' num2str(T)...
        char(10) 'Dyadic negativity: ' num2str(Lambda)...
        char(10) 'Noiseless rotation extent: ' num2str(xi_theta)...
        char(10) 'Channel extent upper bound: ' num2str(xi) ...
        char(10) 'Final Q: ' num2str(Q) '.']);
else

    disp(['Iterations: ' num2str(T)...
        char(10) 'Dyadic negativity: ' num2str(Lambda)...
        char(10) 'Noiseless rotation extent: ' num2str(xi_theta)...
       char(10) 'No valid decomposition found. Final Q: ' num2str(Q)]);
end

end

