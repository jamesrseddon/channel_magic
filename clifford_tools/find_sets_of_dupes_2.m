%%%%% Script to check cleaned Clifford list doesn't have the duplicates
%%%%% any more;

clean_set = two_qubit_list_deduped;
input_set = two_qubit_cliffords;

known_duplicates = two_qubit_duplicates_2copies;

num_to_check = size(two_qubit_duplicates_2copies,1);

dupe_sets = {}
for mm = 1:num_to_check
    this_unitary = input_set{known_duplicates(mm),2};
    [result_vector,matched_rows] = duplicateCheck(this_unitary,clean_set);
    dupe_sets = [dupe_sets;matched_rows];
end
    