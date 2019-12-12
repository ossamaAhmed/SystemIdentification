function [p1_theta_est1,p1_Phi,p1_theta_est2,p1_y_pred] = HS2019_SysID_final_p1_18936872()
%% Solution for Problem 1
clc;
clear all;
close all;
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);
figure_num = 1;
%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data
% Extract Legi from Filename
name=mfilename;
LegiNumber= str2num(name(end-7:end));

[p1_u, p1_y, p1_theta_hat, p1_u_past, p1_y_past,p1_pred_err] = HS2019_SysID_final_p1_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variable p1_U to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Task 1: Obtain initial estimate

%-----------------------------------------------------------------------
%inspecting the output and input first
%-----------------------------------------------------------------------
figure(figure_num);
figure_num = figure_num + 1;
fig = gcf;
fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
plot(p1_u);
hold on;
plot(p1_y);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Input vs. Output';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$k$ $[samples]$';
axes.XLabel.FontSize = 14;
lgd = legend('$u(k)$', '$y(k)$');
lgd.Interpreter = 'latex';
lgd.FontSize = 12;
lgd.Location = 'southwest';

%-----------------------------------------------------------------------
%checking if the input is periodic and identifying the period if its so
%-----------------------------------------------------------------------
N = length(p1_u);
lags = -(N - 1) : 1 : N - 1;
Ru = autocorrelation_periodic(p1_u.', lags);
figure(figure_num);
figure_num = figure_num + 1;
fig = gcf;
fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
stem(lags, Ru);
axes = gca;
axis tight;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Autocorrelation of Input';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\tau$ $[lags]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$R_{u}(\tau)$';
axes.YLabel.FontSize = 14;
%-----------------------------------------------------------------------
%checking the persistency of excitation
%-----------------------------------------------------------------------

%lets count the non zero frequencies in the input spectrum
p1_U = fft(p1_u.');
phi_U = 1 / N * abs(p1_U) .^ 2;

W = 0 : 2*pi/N : pi;
m = size(W, 2);

figure(figure_num);
figure_num = figure_num + 1;
fig = gcf;
fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
stem(W, phi_U((1:m)));
axis([-inf inf 0 inf]);
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Spectrum of Input';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$omega$ $[\frac{rad}{sample}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$\phi_{u}$';
axes.YLabel.FontSize = 14;
%The input seems to be a prbs signal since it has a flat spectrum which
%means its persistently exciting of order 255

%-----------------------------------------------------------------------
%Getting the BLUE estimate
%-----------------------------------------------------------------------
c1 = 0;
c2 = 0;
z = tf('z');
C_z = 1 + c1 / z + c2 / (z^2);

fprintf('A linear model is derived for BLUE estimate: Y = Phi*Theta + Epsilon\n');
fprintf('where, \n');
fprintf('Y(k) = y(k)\n');
fprintf('Theta = [a1; a2; a3; b1; b2; b3]\n');
fprintf('Phi(k) = [-y(k-1) -y(k-2) -y(k-3) u(k-1) u(k-2) u(k-3)]\n');
fprintf('Epsilon(k) = e(k)\n');
fprintf('Error covariance matrix R is expected value of (Epsilon * Epsilon^T) and is taken from the question formulation:\n');
fprintf('Since c1 and c2 are zeros in our case\n');
fprintf('R(i,j) = 0.9 if i = j\n');
fprintf('R(i,j) = 0.2 if 0 < abs(i - j) <= 3\n');
fprintf('R(i,j) = 0 otherwise\n');
%-----------------------------------------------------------------------
%Getting the Covariance Matrix
%-----------------------------------------------------------------------

%forming covariance matrix R now
N2 = N - 3;
R = zeros(N2);
for i = 1 : N2
    for j = 1 : N2
        if i == j
            R(i,j) = 0.9;
        elseif abs(i - j) <=3
            R(i,j) = 0.2;
        end
    end
end

%now lets confirm if R is symmetric
R_difference = abs(R - R.');
R_difference_max = max(R_difference(:));
R_difference_min = min(R_difference(:));
fprintf('It is confirmed that R is symmetric for checking R - R^T maximum: %.4f and minimum: %.4f\n', R_difference_max, R_difference_min);
min_eigen_value_R = min(eig(R));
fprintf('It is confirmed that R is positive definite for checking the smallest eigen value of R > 0: %.4f\n', min_eigen_value_R);

%-----------------------------------------------------------------------
%Getting the Phi Matrix
%-----------------------------------------------------------------------
%now lets compute the BLUE estimator from slide 9.32

%start by forming phi now
n_parameters = 6;
%cutting off the first 3 elements
Phi = zeros(N2, n_parameters);
Y = p1_y(4:end);
for k = 4 : N
    Phi(k-3,:) = [-p1_y(k-1) -p1_y(k-2) -p1_y(k-3) p1_u(k-1) p1_u(k-2) p1_u(k-3)];
end

%-----------------------------------------------------------------------
%Getting the theta estimate now
%-----------------------------------------------------------------------
R_inverse_Phi = R \ Phi;
Z = R_inverse_Phi / (Phi.' * R_inverse_Phi);
theta_estimated = Z.' * Y;
p1_a_estimated = theta_estimated(1:3)
p1_b_estimated = theta_estimated(4:6)
%-----------------------------------------------------------------------
%Self Validation Now
%-----------------------------------------------------------------------
y_estimated = zeros(N,1);
y_estimated(4:end) = Phi * theta_estimated;
estimation_error = p1_y - y_estimated;
figure(figure_num);
figure_num = figure_num + 1;
fig = gcf;
fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
subplot(2,1,1);
plot(p1_y);
hold on;
plot(y_estimated);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Actual Output vs. Model Output';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$k$ $[samples]$';
axes.XLabel.FontSize = 14;
lgd = legend('$y(k)$', '$y_{\hat{\theta}}(k)$');
lgd.Interpreter = 'latex';
lgd.FontSize = 12;
lgd.Location = 'southwest';
subplot(2,1,2);
plot(estimation_error);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Model Output Error';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$k$ $[samples]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$y(k) - y_{\hat{\theta}}(k)$';
axes.YLabel.FontSize = 14;

estimation_error_var = var(estimation_error(4:end), 1);
fprintf('Estimated Model ouput error variance is: %.4f.\n', estimation_error_var);
fprintf('Estimated Model ouput error is: %.4f.\n', sum(abs(estimation_error(4:end))));


p1_theta_est1 = theta_estimated;
p1_Phi = Phi;

%% Task 2: Improve your estimate
c1 = 1;
c2 = pi/4;
z = tf('z');
C_z = 1 + c1 / z + c2 / (z^2);

fprintf('A linear model is derived for BLUE estimate: Y = Phi*Theta + Epsilon\n');
fprintf('where, \n');
fprintf('Y(k) = y(k)\n');
fprintf('Theta = [a1; a2; a3; b1; b2; b3]\n');
fprintf('Phi(k) = [-y(k-1) -y(k-2) -y(k-3) u(k-1) u(k-2) u(k-3)]\n');
fprintf('Epsilon(k) = e(k)\n');
fprintf('Error covariance matrix R is expected value of (Epsilon * Epsilon^T) and is taken from the question formulation:\n');
fprintf('EXPLAIN HERE\n');

%-----------------------------------------------------------------------
%Getting the Covariance Matrix
%-----------------------------------------------------------------------

%forming covariance matrix R now
N2 = N - 5;
R = zeros(N2);
for i = 1 : N2
    for j = 1 : N2
        if i == j
%             R(i,j) = 0.2*(c1 + c2 + c1 + c1*c2 + c2 + c1*c2) + 0.9*(1 + c2^2 + c1^2);
            R(i,j) = 0.9*(1 + c2^2 + c1^2) + 0.2*(c1 + c2 + c1 + c1*c2 + c2 + c1*c2);
        elseif abs(i - j) == 1
%             R(i,j) = 0.2*(1 + c1 + c2 + c2^2 + c1^2 + c1*c2 + c2) + 0.9*(c1 + c1*c2);
            R(i,j) = 0.2*(1 + c1 + c2 + c2^2 + c1^2 + c1*c2 + c2) + 0.9*(c1 + c1*c2);
        elseif abs(i - j) == 2
%             R(i,j) = 0.2*(1 + c1 + c1 + c1^2 + c1*c2 + c1*c2 + c2^2) + c2*0.9;
            R(i,j) = 0.2*(1 + c1 + c1 + c1^2 + c1*c2 + c1*c2 + c2^2) + c2*0.9;
        elseif abs(i - j) == 3
            R(i,j) = (1 + c1 + c1^2 + c2 + c1*c2 + c2^2) * 0.2;
        elseif abs(i - j) == 4
            R(i,j) = (c1 + c2 + c1*c2) * 0.2;
        elseif abs(i - j) == 5
            R(i,j) = c2 * 0.2;
        end
    end
end

%now lets confirm if R is symmetric
R_difference = abs(R - R.');
R_difference_max = max(R_difference(:));
R_difference_min = min(R_difference(:));
fprintf('It is confirmed that R is symmetric for checking R - R^T maximum: %.4f and minimum: %.4f\n', R_difference_max, R_difference_min);
min_eigen_value_R = min(eig(R));
fprintf('It is confirmed that R is positive definite for checking the smallest eigen value of R > 0: %.4f\n', min_eigen_value_R);

%-----------------------------------------------------------------------
%Getting the Phi Matrix
%-----------------------------------------------------------------------
%now lets compute the BLUE estimator from slide 9.32

%start by forming phi now
n_parameters = 6;
%cutting off the first 3 elements
Phi = zeros(N2, n_parameters);
Y = p1_y(6:end);
for k = 6 : N
    Phi(k-5,:) = [-p1_y(k-1) -p1_y(k-2) -p1_y(k-3) p1_u(k-1) p1_u(k-2) p1_u(k-3)];
end

%-----------------------------------------------------------------------
%Getting the theta estimate now
%-----------------------------------------------------------------------
R_inverse_Phi = R \ Phi;
Z = R_inverse_Phi / (Phi.' * R_inverse_Phi);
theta_estimated = Z.' * Y;
p1_a_estimated = theta_estimated(1:3)
p1_b_estimated = theta_estimated(4:6)
%-----------------------------------------------------------------------
%Self Validation Now
%-----------------------------------------------------------------------
y_estimated = zeros(N,1);
y_estimated(6:end) = Phi * theta_estimated;
estimation_error = p1_y - y_estimated;
figure(figure_num);
figure_num = figure_num + 1;
fig = gcf;
fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
subplot(2,1,1);
plot(p1_y);
hold on;
plot(y_estimated);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Actual Output vs. Improved Model Output';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$k$ $[samples]$';
axes.XLabel.FontSize = 14;
lgd = legend('$y(k)$', '$y_{\hat{\theta}}(k)$');
lgd.Interpreter = 'latex';
lgd.FontSize = 12;
lgd.Location = 'southwest';
subplot(2,1,2);
plot(estimation_error);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Improved Model Output Error';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$k$ $[samples]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$y(k) - y_{\hat{\theta}}(k)$';
axes.YLabel.FontSize = 14;

estimation_error_var = var(estimation_error(6:end), 1);
fprintf('Improved Estimated Model ouput error variance is: %.4f.\n', estimation_error_var);
fprintf('Improved Estimated Model ouput error is: %.4f.\n', sum(abs(estimation_error(6:end))));



%-----------------------------------------------------------------------
%Getting the theta estimate now
%-----------------------------------------------------------------------
c1 = 1;
c2 = pi/4;
C_z = 1 + c1 / z + c2 / (z^2);
C_z_inverse = 1 / C_z;
C_z_poles = pole(C_z_inverse);
fprintf('C(z)^(-1) might be stable with only a pole at 0 which might not be stable\n');
%now lets try to see if the C_inv is stable
u_f = lsim(C_z_inverse, p1_u);
y_f = lsim(C_z_inverse, p1_y);

%-----------------------------------------------------------------------
%Self Validation Now
%-----------------------------------------------------------------------
phi_f = zeros(N, n_parameters);
phi_f(2,:) = [-y_f(1) 0 0 u_f(1) 0 0];
phi_f(3,:) = [-y_f(2) -y_f(1) 0 u_f(2) u_f(1) 0];
for k = 4 : N
    phi_f(k,:) = [-y_f(k-1) -y_f(k-2)  -y_f(k-3) u_f(k-1) u_f(k-2) u_f(k-3)];
end
theta_estimated_f = phi_f \ y_f;
y_f_estimated = phi_f * theta_estimated_f;
y_estimated = lsim(C_z, y_f_estimated);
estimation_error = p1_y - y_estimated;
estimation_error_var = var(estimation_error(6:end), 1);
fprintf('Improved Estimated Model ouput error variance is: %.4f.\n', estimation_error_var);
fprintf('Improved Estimated Model ouput error is: %.4f.\n', sum(abs(estimation_error(6:end))));

%% Task 3: Compute prediction

p1_y_pred = [];
end
