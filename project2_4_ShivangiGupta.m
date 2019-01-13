clc; clear all; format compact;

%%Problem 4 (Project 2)
%   SHIVANGI GUPTA
clc; clear all; format compact;

load('project2_data.mat');
data = data;

samples = data;
%----------------------ONE COMPONENT MIXTURE-------------------------------
est_sigma1 = 10 * rand(1,1);
est_p1 = 10 * rand(1,1);
est_p1      = est_p1 / sum(est_p1);
l1 = [];

for ii = 1:100

        for jj=1:numel(est_sigma1)
            w(jj,:) = est_p1(jj) * raylpdf(samples, est_sigma1(jj));
        end

        for jj=1:numel(est_sigma1)
             rayl(jj,:) = w(jj,:) ./ sum(w,1);
             N(jj) = sum(rayl(jj,:),2);
             est_sigma1_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
             est_sigma1(jj) = sqrt(est_sigma1_sq(jj));
        end
        
        l1(end+1) = sum(log(sum(w(:,:))));
        
        for jj=1:numel(est_sigma1)
            est_p1(jj) = N(jj)./sum(N(:,:));
        end

        if numel(l1) > 2
            if abs(l1(end) - l1(end - 1)) < 1E-3
                break
            end
        end
end

%----------------------TWO COMPONENT MIXTURE-------------------------------
est_sigma2 = 10 * rand(1,2);
est_p2 = 10 * rand(1,2);
est_p2      = est_p2 / sum(est_p2);
l2 = [];

for ii = 1:100

        for jj=1:numel(est_sigma2)
            w(jj,:) = est_p2(jj) * raylpdf(samples, est_sigma2(jj));
        end

        for jj=1:numel(est_sigma2)
             rayl(jj,:) = w(jj,:) ./ sum(w,1);
             N(jj) = sum(rayl(jj,:),2);
             est_sigma2_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
             est_sigma2(jj) = sqrt(est_sigma2_sq(jj));
        end
        
        l2(end+1) = sum(log(sum(w(:,:))));
        
        for jj=1:numel(est_sigma2)
            est_p2(jj) = N(jj)./sum(N(:,:));
        end

        if numel(l2) > 2
            if abs(l2(end) - l2(end - 1)) < 1E-3
                break
            end
        end
end

%----------------------THREE COMPONENT MIXTURE-------------------------------
est_sigma3 = 10 * rand(1,3);
est_p3 = 10 * rand(1,3);
est_p3      = est_p3 / sum(est_p3);
l3 = [];

for ii = 1:100

        for jj=1:numel(est_sigma3)
            w(jj,:) = est_p3(jj) * raylpdf(samples, est_sigma3(jj));
        end

        for jj=1:numel(est_sigma3)
             rayl(jj,:) = w(jj,:) ./ sum(w,1);
             N(jj) = sum(rayl(jj,:),2);
             est_sigma3_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
             est_sigma3(jj) = sqrt(est_sigma3_sq(jj));
        end
        
        l3(end+1) = sum(log(sum(w(:,:))));
        
        for jj=1:numel(est_sigma3)
            est_p3(jj) = N(jj)./sum(N(:,:));
        end

        if numel(l3) > 2
            if abs(l3(end) - l3(end - 1)) < 1E-3
                break
            end
        end
end

%----------------------FOUR COMPONENT MIXTURE-------------------------------
est_sigma4 = 10 * rand(1,4);
est_p4 = 10 * rand(1,4);
est_p4      = est_p4 / sum(est_p4);
l4 = [];

for ii = 1:100

        for jj=1:numel(est_sigma4)
            w(jj,:) = est_p4(jj) * raylpdf(samples, est_sigma4(jj));
        end

        for jj=1:numel(est_sigma4)
             rayl(jj,:) = w(jj,:) ./ sum(w,1);
             N(jj) = sum(rayl(jj,:),2);
             est_sigma4_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
             est_sigma4(jj) = sqrt(est_sigma4_sq(jj));
        end
        
        l4(end+1) = sum(log(sum(w(:,:))));
        
        for jj=1:numel(est_sigma4)
            est_p4(jj) = N(jj)./sum(N(:,:));
        end

        if numel(l4) > 2
            if abs(l4(end) - l4(end - 1)) < 1E-3
                break
            end
        end
end

%----------------------FIVE COMPONENT MIXTURE-------------------------------
est_sigma5 = 10 * rand(1,5);
est_p5 = 10 * rand(1,5);
est_p5      = est_p5 / sum(est_p5);
l5 = [];

for ii = 1:100

        for jj=1:numel(est_sigma5)
            w(jj,:) = est_p5(jj) * raylpdf(samples, est_sigma5(jj));
        end

        for jj=1:numel(est_sigma5)
             rayl(jj,:) = w(jj,:) ./ sum(w,1);
             N(jj) = sum(rayl(jj,:),2);
             est_sigma5_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
             est_sigma5(jj) = sqrt(est_sigma5_sq(jj));
        end
        
        l5(end+1) = sum(log(sum(w(:,:))));
        
        for jj=1:numel(est_sigma5)
            est_p5  (jj) = N(jj)./sum(N(:,:));
        end

        if numel(l5) > 2
            if abs(l5(end) - l5(end - 1)) < 1E-3
                break
            end
        end
end

%----------------------SIX COMPONENT MIXTURE-------------------------------
est_sigma6 = 10 * rand(1,6);
est_p6 = 10 * rand(1,6);
est_p6      = est_p6 / sum(est_p6);
l6 = [];

for ii = 1:100

        for jj=1:numel(est_sigma6)
            w(jj,:) = est_p6(jj) * raylpdf(samples, est_sigma6(jj));
        end

        for jj=1:numel(est_sigma6)
             rayl(jj,:) = w(jj,:) ./ sum(w,1);
             N(jj) = sum(rayl(jj,:),2);
             est_sigma6_sq(jj) = sum(rayl(jj, :) .* (samples.').^2) / (2.*N(jj));
             est_sigma6(jj) = sqrt(est_sigma6_sq(jj));
        end
        
        l6(end+1) = sum(log(sum(w(:,:))));
        
        for jj=1:numel(est_sigma6)
            est_p6(jj) = N(jj)./sum(N(:,:));
        end

        if numel(l6) > 2
            if abs(l6(end) - l6(end - 1)) < 1E-3
                break
            end
        end
end

%-------------------------------HISTOGRAM----------------------------------

pdf1 = @(input) arrayfun(@(x) dot(est_p1, raylpdf(x, est_sigma1)), input);   
pdf2 = @(input) arrayfun(@(x) dot(est_p2, raylpdf(x, est_sigma2)), input);
pdf3 = @(input) arrayfun(@(x) dot(est_p3, raylpdf(x, est_sigma3)), input);
pdf4 = @(input) arrayfun(@(x) dot(est_p4, raylpdf(x, est_sigma4)), input);
pdf5 = @(input) arrayfun(@(x) dot(est_p5, raylpdf(x, est_sigma5)), input);
pdf6 = @(input) arrayfun(@(x) dot(est_p6, raylpdf(x, est_sigma6)), input);

fprintf('One Component mixture Model \n');
fprintf('Estimated value of sigma1 : \n' );
fprintf('%f \n',est_sigma1(:));
fprintf('\nEstimated value of p1 : \n' );
fprintf('%f \n',est_p1(:));

fprintf('Two Component mixture Model \n');
fprintf('Estimated value of sigma2 : \n' );
fprintf('%f \n',est_sigma2(:));
fprintf('\nEstimated value of p2 : \n' );
fprintf('%f \n',est_p2(:));

fprintf('Three Component mixture Model \n');
fprintf('Estimated value of sigma3 : \n' );
fprintf('%f \n',est_sigma3(:));
fprintf('\nEstimated value of p3 : \n' );
fprintf('%f \n',est_p3(:));

fprintf('Four Component mixture Model \n');
fprintf('Estimated value of sigma4 : \n' );
fprintf('%f \n',est_sigma4(:));
fprintf('\nEstimated value of p4 : \n' );
fprintf('%f \n',est_p4(:));

fprintf('Five Component mixture Model \n');
fprintf('Estimated value of sigma5 : \n' );
fprintf('%f \n',est_sigma5(:));
fprintf('\nEstimated value of p5 : \n' );
fprintf('%f \n',est_p5(:));

fprintf('Six Component mixture Model \n');
fprintf('Estimated value of sigma6 : \n' );
fprintf('%f \n',est_sigma6(:));
fprintf('\nEstimated value of p6 : \n' );
fprintf('%f \n',est_p6(:));

figure(1)
histogram(samples, 'Normalization', 'pdf');
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf1(t), 'LineWidth', 2)
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf2(t), 'LineWidth', 2)
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf3(t), 'LineWidth', 2)
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf4(t), 'LineWidth', 2)
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf5(t), 'LineWidth', 2)
hold on
t = linspace(min(samples), max(samples), 1000);
plot(t, pdf6(t), 'LineWidth', 2)
hold off

legend('Sample Data','One component mixture', 'Two component mixture', 'Three component mixture',...
    'Four component mixture', 'Five component mixture', 'Six component mixture')
   
