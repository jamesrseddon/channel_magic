function [] = syntheticNoiseAnalysis(result_set,T,delta,N,T_N)
%SYNTHETICNOISEANALYSIS Analyse results from synthetic noise sweep.
%   (1) Plot mixing parameter lambda against dyadic negativity.
%   (2) Plot total error for some fixed number of samples T and failure
%   probability delta.
lambda_vector = cell2mat(result_set(2:end,1));
dyad_neg_vector = cell2mat(result_set(2:end,3));

lambda_label = result_set{1,1};
dyadic_label = result_set{1,3};

figure
plot(lambda_vector,dyad_neg_vector,'o-');
xlabel(lambda_label);
ylabel(dyadic_label);


product_vec = lambda_vector.*dyad_neg_vector;

figure
plot(lambda_vector,product_vec,'o-');
xlabel(lambda_label);
ylabel(dyadic_label);
title('product of lambda and dyadic negativity');


[result_rows,~] = size(result_set);

M = result_rows - 1;

const = sqrt(log(1/delta)/T);

Delta = zeros(1,M);

plot_title = ['number of samples = '...
                num2str(T) ', success probability = ' ...
                num2str(1 - delta)];

for kk = 1:M
    lambda = lambda_vector(kk);
    dyad_neg = dyad_neg_vector(kk);
    Delta(kk) = (lambda - 1) + lambda*dyad_neg*const;
end

figure
plot(lambda_vector,Delta,'o-g');
xlabel(lambda_label);
ylabel('Total error bound');
title(plot_title);

% many copies

Delta_tot = zeros(1,M);
lambda_tot = lambda_vector.^N;
dyad_neg_tot = dyad_neg_vector.^N;
product_vec = lambda_tot.*dyad_neg_tot;
const_N = sqrt(log(1/delta)/T_N);
dyad_neg_max = dyad_neg_tot(1);
const_norm = 1/(2*dyad_neg_max);


figure
semilogy(lambda_vector,product_vec,'o-');
xlabel(lambda_label);
ylabel(dyadic_label);
title('product of lambda^N and dyadic negativity^N');


for kk = 1:M
    lambda = lambda_tot(kk);
    dyad_neg = dyad_neg_tot(kk);
    Delta_tot(kk) = (lambda - 1) + lambda*dyad_neg*const_N;
    Delta_tot_norm(kk) = (lambda - 1) + lambda*dyad_neg*const_norm;
end

figure
semilogy(lambda_vector,Delta_tot,'o-g');
xlabel(lambda_label);
ylabel('Total error bound');
title(['Many copies, N = ' num2str(N) ', number of samples = ' num2str(T_N)]);

figure
semilogy(lambda_vector,Delta_tot_norm,'o-g');
xlabel(lambda_label);
ylabel('Total error bound');
title(['Many copies normalised, N = ' num2str(N)]);




