function [member_flag] = checkMember(A,S,tolerance)
%CHECKMEMBER Check if matrix A is in set S.
%   S should be specified as 3-dimensional array where third dimension
%   indexes the element of the set, and dimensions 1 and 2 are the
%   dimensions of a square matrix.
%   ouput is True if A is in S, up to tolerance, False otherwise.
member_flag = false;

if isempty(S)
    return
end

num_elements = size(S,3);


for kk = 1:num_elements
    B = S(:,:,kk);
    if checkEqual(A,B,tolerance)
        member_flag = true;
        return
    end
end

