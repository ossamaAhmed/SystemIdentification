% Initialization
clc;
clear all;
close all;

% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%get the data
load('Data_ex9.mat');
y_true = ex9_y;
u_original = ex9_u;

num_of_parameters = 3;
N = length(y_true);
n_parameters = 4;
fun = @(x)fminconObjective(x, n_parameters);
x0 = zeros(n_parameters + N, 1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @(x)fminconConstraint(x, y_true, u_original);
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000);
x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
a1 = x(1)
a2 = x(2)
b1 = x(3)
c1 = x(4)
%get the error
e = x(5:end);
y_est = y_true - e;
figure(1);
fig = gcf;
subplot(2,1,1);
plot(y_true);
hold on;
plot(y_est);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'y Actual vs. Estimate';
axes.Title.FontSize = 18;
% Error
subplot(2,1,2);
plot(y_true - y_est);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Self Validation Error';
axes.Title.FontSize = 18;

%-------------------------------------------------------
%-------------------PART TWO---------------------------
%-------------------------------------------------------
phi = zeros(N, 3);
phi(2,:) = [-y_true(1) 0 u_original(1)];
for k = 3 : N
    phi(k, :) = [-y_true(k-1) -y_true(k-2) u_original(k-1)];
end

%using LS estimate
thetaLS = phi \ y_true;
a1LS = thetaLS(1)
a2LS = thetaLS(2)
b1LS = thetaLS(3)
z = tf('z');
G = (b1LS / z) / (1 + a1LS / z + a2LS / z^2);
x = lsim(G, u_original);

eta = zeros(3, N);
Rk = zeros(3);
eta(:,2) = [-x(1); 0; u_original(1)];
Rk = Rk + eta(:,2) * phi(2,:);

eta_y = zeros(3, 1);
eta_y = eta_y + eta(:,2) * y_true(2);
for k = 3 : N
    eta(:,k) = [-x(k-1); -x(k-2); u_original(k-1)];
    Rk = Rk + eta(:,k) * phi(k,:);
    eta_y = eta_y + eta(:,k) * y_true(k);
end

Rk = Rk / N;
eta_y = eta_y / N;

thetaIV = Rk \ eta_y;
a1IV = thetaIV(1)
a2IV = thetaIV(2)
b1IV = thetaIV(3)

