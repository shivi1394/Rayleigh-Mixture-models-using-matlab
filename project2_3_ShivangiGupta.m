clc; clear all; format compact;

%%Problem 3 (Project 2)
%   SHIVANGI GUPTA
clc; clear all; format compact;

sigma = [2 5];
p = [2 1];
p = p / sum(p);


latent_variable = sum(rand(1E5, 1) > cumsum(p), 2) + 1;

pdf = @(input) arrayfun(@(x) dot(p, raylpdf(x, sigma)), input);

group1_samples = raylrnd(sigma(1), sum(latent_variable == 1), 1);
group2_samples = raylrnd(sigma(2), sum(latent_variable == 2), 1);
samples  = [group1_samples; group2_samples];

%--------------------------------------------------------------------------
%EXPECTATION MAXIMIZATION
est_sigma = 10 * rand(1,2);

est_p = 10 * rand(1,2);
est_p = est_p / sum(est_p);
l = [];

for ii = 1:100
    
    for jj=1:numel(est_sigma)
        w(jj,:) = est_p(jj) * raylpdf(samples, est_sigma(jj));
    end
    
    for jj=1:numel(est_sigma)
         rayl(jj,:) = w(jj,:) ./ sum(w,1);
         N(jj) = sum(rayl(jj,:),2);
         est_sigma_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
         est_sigma(jj) = sqrt(est_sigma_sq(jj));
    end
    
    l(end+1) = sum(log(sum(w(:,:))));
    
    for jj=1:numel(est_sigma)
        est_p(jj) = N(jj)./sum(N(:,:));
    end
    
    if numel(l) > 2
        if abs(l(end) - l(end - 1)) < 1E-3
            break
        end
    end
end


pdf1 = @(input) arrayfun(@(x) dot(est_p, raylpdf(x, est_sigma)), input);
figure(1)
histogram(samples, 'Normalization', 'pdf')
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf(t),'r', 'LineWidth', 2)
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf1(t),'b', 'LineWidth', 2)
hold off

legend('Sample Data','pdf of prescribed mixture','Arbitrary mixture model');

fprintf('Estimated value of sigma : \n' );
fprintf('%f \n',est_sigma(:));
fprintf('\nEstimated value of p : \n' );
fprintf('%f \n',est_p(:));

