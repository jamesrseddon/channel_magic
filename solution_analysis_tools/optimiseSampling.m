function [sampling_opt,lambda_opt] = optimiseSampling(Delta,result_set)
%OPTIMISESAMPLING Plot required sampling constant to achieve given error in
% synthetic noise algorithm.
%   result_table is the result from a synthetic noise sweep.
%   Delta is the desired error bound.

lambda_vector = cell2mat(result_set(2:end,1));
cost_vector = cell2mat(result_set(2:end,3));


lambda_label = result_set{1,1};
cost_label = result_set{1,3};

num_points = length(lambda_vector)

K = zeros(num_points,1);

for kk = 1:num_points
    lambda = lambda_vector(kk);
    negativity = cost_vector(kk);
    K(kk) = (Delta + 1 - lambda)/(negativity*lambda);
end

figure
plot(lambda_vector,K,'.-')
xlabel(lambda_label)
ylabel('sampling constant')
title(['Required sampling constant for max error = ' num2str(Delta) '.'])

[sampling_opt opt_index] = max(K);

lambda_opt = lambda_vector(opt_index);

end

