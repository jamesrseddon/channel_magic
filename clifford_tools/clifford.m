function output = clifford(gatestring)
%%% build the single qubit clifford gates S and H, or the CNOT.
%%% outputs the unitary in matrix form, based on input string:
%%% gatestring = 'S','H' or 'CNOT'

cliffordMap = containers.Map();

cliffordMap('S')=[ 1 0 ;
                   0 1i];
cliffordMap('H')= 1/(sqrt(2))*[ 1 1;
                        1 -1];
                    
cliffordMap('CNOT') = [1 0 0 0;
                       0 1 0 0;
                       0 0 0 1;
                       0 0 1 0];
                   
output = cliffordMap(gatestring);