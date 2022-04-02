function [ clifford_set ] = generateCliffordRotations()
%GENERATECLIFFORDROTATIONS Generate the SO(3) representations of the
%single-qubit Clifford gates.


R_S = [0 -1 0;
       1 0 0;
       0 0 1];

R_H = [0 0 1;
       0 -1 0;
       1 0 0];

current_set = zeros(3,3,2);

current_set(:,:,1) = R_S;
current_set(:,:,2) = R_H;

incomplete = true;

while incomplete
    N = size(current_set,3);
    num_new_elements = 0;
for jj = 1:N
    for kk = 1:N
        candidate = current_set(:,:,jj)*current_set(:,:,kk);
        found_match = false;
        for mm = 1:(N + num_new_elements)
            if checkEqual(candidate,current_set(:,:,mm),1e-12)
                found_match = true;
                break
            end
        end
        if ~found_match
            num_new_elements = num_new_elements + 1;
            new_index = N + num_new_elements;
            current_set(:,:,new_index) = candidate;
        end
    end
end
    disp(['Number new elements: ' num2str(num_new_elements)]);
    if num_new_elements == 0
        incomplete = false;
    end
end
total_elements = size(current_set,3);

disp(['Total elements: ' num2str(total_elements)]);
clifford_set = current_set;