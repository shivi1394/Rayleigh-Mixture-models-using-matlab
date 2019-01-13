clc; clear all; format compact;

%%Problem 1 (Project 2)
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

figure(1)
histogram(samples, 'Normalization', 'pdf')
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf(t), 'LineWidth', 2)
hold off
grid on

