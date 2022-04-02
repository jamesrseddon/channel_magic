function [equivalence_vector,matchrows] = duplicateCheck(test_matrix,test_set)
test_matrix;

test_set;

set_size = size(test_set,1);

counter = 0;
matchrows = [];
for kk=1:set_size
    equivalent(kk) = phaseEquiv(test_matrix,test_set{kk,2});
    if equivalent(kk)
        display(['Found a match at row ' num2str(kk) '.']);
        matchrows = [matchrows kk];
        counter = counter + 1;
        display(['Count so far = ' num2str(counter)]);
    end
    if rem(kk,100) == 0
        % display(['Now at row ' num2str(kk)])
    end
end

display(['Final count = ' num2str(sum(equivalent))]);

equivalence_vector = equivalent;