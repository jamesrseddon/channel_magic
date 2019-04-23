For the files 3QubitCalculations.txt and 4QubitCalculations.txt the column headings 
{"Ind", "T_cost_AM", "R", "circ"}
correspond to 
Ind: Index
T_cost_AM: The T cost of synthesizing the unitary over the CNOT+T gate set using arXiv:1601.07363 [quant-ph]
R: Robustness of the state |U>
circ: Description of the circuit implementing U in terms of T gates, CS gates and CCZ gates. Numbers after T or CS/CCZ denote the qubits or pairs/triples of qubits respectively, on which the gate acts. For example "T{1, 2, 3}CS{1, 3, 2, 3}CCZ{}" means apply T on qubits 1,2 and 3 and CS on qubits "1 and 3" and "2 and 3".



The files Amat2.mat ,..., Amat5.mat correspond to the A matrices for the linear system Ax=b. Each column corresponds to a stabilizer state written in the Pauli operator basis. The same ordering should be used for decomposing a state rho as vector b. We use the ordering
IIIII
IIIIX
IIIIY
IIIIZ
IIIXI
IIIXX
.....
ZZZZZ

The syntax for calculating robustness (assuming a working CVX [http://cvxr.com/cvx/] installation) is :

%n is one of {60,1080,36720,2423520}
cvx_begin 
    variable x(n)
    minimize( norm( x, 1 ) )
    subject to
       A * x == b
cvx_end
robustness=cvx_optval;

% another option is to ask for a certificate, y, of optimality
cvx_begin
    variable x(n)
    dual variable y;
    minimize( norm( x, 1 ) )
    subject to
       y: A * x == b
cvx_end
robustness=cvx_optval;
