% Initialization
clc;
clear all;
close all;

% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%define the model now
z = tf('z');
A_z = 1 - (1.5 / z) + (0.7 / (z^2));
B_z = (1.0 / z) + (0.5 / (z^2));
C_z = 1 - (1 / z) + (0.2 / (z^2));
part_1 = 0;
%generate e_u
var = 1;
std_dev = sqrt(var);

%define the input now
N = 10^4;
e_u = std_dev * randn(N, 1);
u = zeros(N, 1);
u(2) = 0.1 * u(1) + e_u(1);
for k = 3 : N
    u(k) = 0.1 * u(k-1) + 0.12 * u(k-2) + e_u(k-1) + 0.2 * e_u(k-2);
end

%obtain now least squares estimate
C_z_inv = 1 / C_z;
u_f = lsim(C_z_inv, u);

%define the parameter vector theta
a1 = -1.5;
a2 = 0.7;
b1 = 1;
b2 = 0.5;
c1 = -1;
c2 = 0.2;
C_params = [c2 c1 1];
theta = [a1; a2; b1; b2];
num_of_experiments = 60;
a_1_est = zeros(1, num_of_experiments);
a_2_est = zeros(1, num_of_experiments);
b_1_est = zeros(1, num_of_experiments);
b_2_est = zeros(1, num_of_experiments);
num_of_parameters = size(theta, 1);
y_true = zeros(N, 1);
y_est = zeros(N, 1);
e_filtered = zeros(N);
phi = zeros(N, num_of_parameters);
phi_f = zeros(N, num_of_parameters);
for i = 1 : num_of_experiments
    %generating noise
    e = std_dev * randn(N, 1);
    %generate y_true
    u_transformed = lsim(B_z/A_z, u);
    e_transformed = lsim(C_z/A_z, e);
    y_true = u_transformed + e_transformed;
    %estimate the params now
    y_f = lsim(C_z_inv, y_true);
    phi_f(2,:) = [-y_f(1) 0 u_f(1) 0];
    for k = 3 : N
        phi_f(k,:) = [-y_f(k-1) -y_f(k-2) u_f(k-1) u_f(k-2)];
    end
    theta_est = phi_f \ y_f;
    a_1_est(i) = theta_est(1);
    a_2_est(i) = theta_est(2);
    b_1_est(i) = theta_est(3);
    b_2_est(i) = theta_est(4);
    
    %create phi
    y_f_est = phi_f * theta_est;
    y_est = lsim(C_z, y_f_est);
    if part_1
        figure(1);
        fig = gcf;
        % y Actual vs. Estimate
        subplot(2,1,1);
        plot(y_true);
        hold on;
        plot(y_est);
        axis tight;
        axes = gca;
        axes.Title.Interpreter = 'latex';
        axes.Title.String = 'Self Validation y Actual vs. Estimate';
        axes.Title.FontSize = 18;
        % Error
        subplot(2,1,2);
        plot(y_true - y_est);
        axis tight;
        axes = gca;
        axes.Title.Interpreter = 'latex';
        axes.Title.String = 'Self Validation Error';
        axes.Title.FontSize = 18;
    end
end
a_1_est_avg = mean(a_1_est);
a_2_est_avg = mean(a_2_est);
b_1_est_avg = mean(b_1_est);
b_2_est_avg = mean(b_2_est);
figure(2);
subplot(2,2,1);
histogram(a_1_est);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'a1 Estimate';
axes.Title.FontSize = 18;
% a2 Histogram
subplot(2,2,2);
histogram(a_2_est);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'a2 Estimate';
axes.Title.FontSize = 18;
% b1 Histogram
subplot(2,2,3);
histogram(b_1_est);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'b1 Estimate';
axes.Title.FontSize = 18;
% a2 Histogram
subplot(2,2,4);
histogram(b_2_est);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'b2 Estimate';
axes.Title.FontSize = 18;


