function [expected_overlap,overlap_list] = expectedOverlap(distribution,STAB_mat)
%EXPECTEDOVERLAP  Evaluate <psi|E(|w><w|)|psi> where |w> are terms in the
%stabiliser decomposition of |psi>.
%   We take a pure state |psi>, decomposed into stabiliser states as per
%   BBCCGH:
%      |psi> = SUM_j c_j |phi_j>.
%   The state |w> is a random variable where |w> = |W_j> with probability
%   p_j = |c_j|/||c||_1, and where |W_j> = (c_j/|c_j|)|phi_j>.


[rows,cols] = size(STAB_mat);

abs_dist = abs(distribution);
one_norm = sum(abs_dist);
prob_dist = abs_dist/one_norm;
phase_dist = distribution./abs_dist;

psi = STAB_mat*distribution;

overlap_list = abs(STAB_mat'*psi).^2;

overlap_prob = prob_dist.*overlap_list;

expected_overlap = sum(overlap_prob);

end

