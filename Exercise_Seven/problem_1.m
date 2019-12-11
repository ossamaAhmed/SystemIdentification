% Initialization
clc;
clear all;
close all;

% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%define the true model
a = 1/2;
b = 1;
p = 2;

%define the random signal proporties
var = 0.2;
std_dev = sqrt(var);
low = -sqrt(12) / 2 * std_dev;
upper = sqrt(12) / 2 * std_dev;
orders = 13;
num_of_experiments = 1000;
a_normal_est = zeros(1, num_of_experiments);
b_normal_est = zeros(1, num_of_experiments);
a_uniform_est = zeros(1, num_of_experiments);
b_uniform_est = zeros(1, num_of_experiments);
SSR = zeros(1, num_of_experiments);
figure_num = 1;
figure(1);
for i = 5 : orders
    N = 2 ^ (i + 1);
    %fix an input
    u = randn(N, 1);
    y_normal_true = zeros(N, 1);
    y_uniform_true = zeros(N, 1);
    phi_normal = zeros(N, 2);
    phi_uniform = zeros(N, 2);
    for j = 1 : num_of_experiments
        w_normal = std_dev * randn(N, 1);
        w_uniform = low + (upper - low) * rand(N, 1);
        y_normal_true(1) = w_normal(1);
        y_uniform_true(1) = w_uniform(1);
        for k = 2 : N
            y_normal_true(k) = a * y_normal_true(k-1) + b * u(k-1) + w_normal(k);
            phi_normal(k,:) = [-y_normal_true(k-1) u(k-1)];
            y_uniform_true(k) = a * y_uniform_true(k-1) + b * u(k-1) + w_uniform(k);
            phi_uniform(k,:) = [-y_uniform_true(k-1) u(k-1)];
        end
        theta_normal = phi_normal \ y_normal_true;
        a_normal_est(j) = -theta_normal(1);
        b_normal_est(j) = theta_normal(2);
        theta_uniform = phi_uniform \ y_uniform_true;
        a_uniform_est(j) = -theta_uniform(1);
        b_uniform_est(j) = theta_uniform(2);
        
        %calculating residuals
        y_estimated = phi_normal * theta_normal;
        SSR(j) = sum((y_normal_true - y_estimated).^2) / var;
    end
    SSRmin = min(SSR);
    SSRmax = max(SSR);
    SSRvec = SSRmin : (SSRmax - SSRmin) / num_of_experiments : SSRmax;
    SSRchi2 = chi2pdf(SSRvec, N - p);
    SSRnorm = normpdf(SSRvec, N - p, sqrt(2 * (N - p)));
    
    subplot(orders,5, ((i-5)*5)+1);
    histogram(a_normal_est);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = strcat('Normal a Estimate for N = ', int2str(N));
    axes.Title.FontSize = 18;
    % Normal b Histogram
    subplot(orders,5, ((i-5)*5)+2);
    histogram(b_normal_est);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = strcat('Normal b Estimate for N = ', int2str(N));
    axes.Title.FontSize = 18;
    % Unifrom a Histogram
    subplot(orders,5, ((i-5)*5)+3);
    histogram(a_uniform_est);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = strcat('Unifrom a Estimate for N = ', int2str(N));
    axes.Title.FontSize = 18;
    % Unifrom b Histogram
    subplot(orders,5, ((i-5)*5)+4);
    histogram(b_uniform_est);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = strcat('Uniform b Estimate for N = ', int2str(N));
    axes.Title.FontSize = 18;
    % Normal SSR Histogram
    subplot(orders,5, ((i-1)*5)+5);
    histogram(SSR, 'Normalization', 'pdf');
    hold on;
    plot(SSRvec, SSRchi2);
    plot(SSRvec, SSRnorm);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = strcat('Normal SSR Estimate for N = ', int2str(N));
    axes.Title.FontSize = 18;

end
