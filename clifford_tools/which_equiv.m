%%% check through a set of unitaries to see if it contains any duplicates
clear occurrence_count duplicate_count

set_to_test = two_qubit_list_deduped;
how_many = 500;

set_size = size(set_to_test,1)
duplicate_count = 0; 
tic
for uu = 1:how_many
    display(['Iteration ' num2str(uu) ' of ' num2str(how_many) '.']);
    unitary_index = set_size + 1 - uu;
    this_unitary = set_to_test{unitary_index,2};
    equivalence_vector = duplicateCheck(this_unitary,set_to_test);
    occurrence_count(uu,1) = unitary_index;
    occurrence_count(uu,2) = sum(equivalence_vector);
    if occurrence_count(uu,2) > 1
        display(['Unitary number ' num2str(unitary_index) ' occurs '...
            num2str(occurrence_count(uu,2)) ' times!']);
        duplicate_count = duplicate_count + 1;
    end
end
display(['Total duplicates spotted : ' num2str(duplicate_count)]); 

toc
    