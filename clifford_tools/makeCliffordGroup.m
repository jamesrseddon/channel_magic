function output = makeCliffordGroup(n);
%%% Take generators and enumerate all the unique unitaries in the Clifford
%%% group for n qubits.

if ~(rem(n,1)==0 && n > 0)
    msg = 'n must be a positive integer';
    error(msg);
end

dim = 2^n;

%%% Define the generators and the identity
S = clifford('S');
H = clifford('H');
I = eye(2);

I_n = eye(dim);

% ****Need to add CNOTs and tensor products for the multiqubit gates****
CNOT = clifford('CNOT');

%%% Calculate number of elements to find

current_count = 1;

for kk=1:n
    current_count = 2*(4^kk -1)*4^kk*current_count;
end

C_size = current_count*10;

%%% Make the set of generators
gens = {};

if n == 1
    gens = {'S', S;
        'H', H;};
end

%%% make the single-qubit generators
if n > 1
    for mm = 1:n
        S_label = ['S(' num2str(mm) ')'];
        H_label = ['H(' num2str(mm) ')'];
        S_matrix = 1;
        H_matrix = 1;
        for kk = 1:n
            if kk == mm
                S_matrix = kron(S_matrix,S);
                H_matrix = kron(H_matrix,H); 
            else
                S_matrix = kron(S_matrix,I);
                H_matrix = kron(H_matrix,I); 
            end
        end
        new_generator_rows = { S_label, S_matrix;
                           H_label, H_matrix;};
        gens = [gens; new_generator_rows];
    end
end
%%% make the CNOT generators
if n > 1
    %%% make a CNOT with control 1 target 2. This will be permuted to get
    %%% other CNOTs.
    base_CNOT = CNOT;
    if n > 2
        dim_rem = 2^(n-2);
        I_rem = eye(dim_rem);
        base_CNOT = kron(CNOT,I_rem);
    end
    
    sub_dims = 2*ones(1,n); % set subsystem dimensions.
    
    %%% loop over controls
    for mm = 1:(n-1)
        %%% loop over targets
        for kk = (mm+1):n
            C_label = ['CNOT(' num2str(mm) '>' num2str(kk) ')'];
            %%% determine which positions we need to switch
            permuter = [1:n];
           
            %%% switch the target qubit index
            %%% this needs to be done first, as otherwise can permute
            %%% incorrectly when the control qubit is 2.
            permuter([2,kk]) = permuter([kk,2]);
            %%% switch the control qubit index
            permuter([1,mm]) = permuter([mm,1]); 
           
            %%% permute CNOT(1>2) to CNOT(mm>kk)
            C_matrix = syspermute(base_CNOT,permuter,sub_dims);
            new_row = {C_label C_matrix};
            gens = [gens; new_row];
        end
    end
    
end




num_gens = size(gens,1);
    
%%% Make the initial set

current_set = {'I', I_n};

current_set = [current_set; gens];
test_set = gens;
set_size = size(current_set,1);
iterations = 0

%%% Loop to produce a new set of elements
while and(set_size < C_size, iterations < 5000)
    size_last_run = set_size;
    iterations = iterations + 1;
    display(['Iteration number ' num2str(iterations)]);
    num_to_check = size(test_set);
    new_test_set = {};
    for pp = 1:num_to_check;
        for qq = 1:2*num_gens;
            if qq<=num_gens
                % left multiply
                new_name = [gens{qq,1} test_set{pp,1}];
                new_matrix = gens{qq,2}*test_set{pp,2};
            else
                % right multiply - probably not necessary?
%                 new_name = [test_set{pp,1} gens{qq-num_gens,1}];
%                 new_matrix = test_set{pp,2}*gens{qq-num_gens,2};
            end
            found_duplicate = false; % flag for checking if already have it
            set_size = size(current_set,1);
            % compare with current list to see if we already have it
            for rr = 1:set_size
                equivalence_check = phaseEquiv(new_matrix,...
                                                    current_set{rr,2});
                if equivalence_check
                    found_duplicate = true;
                    break
                end
            end
            % if we don't already have it, add it to the set.
            if ~found_duplicate
                new_row = { new_name, new_matrix };
                current_set = [current_set; new_row];
                % also add it to the test set for the next run
                new_test_set = [new_test_set; new_row];
            end
        end
    end
    test_set = new_test_set;
    set_size = size(current_set,1);
    display(['Current set size : ' num2str(set_size)]);
    found_this_run = set_size - size_last_run;
    display(['Unitaries added this time: ' num2str(found_this_run)]);
    if found_this_run == 0;
        break
    end
end

output= current_set;