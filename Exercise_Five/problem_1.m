% Initialization
clc;
clear all;
close all;

% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%define the true model
gamma = 5 ;
true_y = generate_prbs(gamma, [0, 1, 0, 1, 0, 1]);
figure(1);
plot(true_y);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('True Signal Y');
axes.Title.FontSize = 18;
%---------------------------------------------------------------------
%------------------------PART TWO-------------------------------------
%---------------------------------------------------------------------
M = 2^6 -1 ; 
y_period = true_y(1:M) ;
lags = 0:M-1;
autocorrelation_y_periodic = autocorrelation_periodic(y_period, lags);
autocorrelation_y_finite = autocorrelation_finite(y_period, lags);
figure(2);
subplot(2, 1, 1);
plot(0:M-1, autocorrelation_y_periodic);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('Autocorrelation periodic');
axes.Title.FontSize = 18;
subplot(2, 1, 2);
plot(0:M-1, autocorrelation_y_finite);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('Autocorrelation finite');
axes.Title.FontSize = 18;
%---------------------------------------------------------------------
%------------------------PART THREE-------------------------------------
%---------------------------------------------------------------------

%generate y(k) with different initial conditions
x_initial = eye(6);
x_initial(7, :) = [1, 1, 0, 0, 1, 1];
x_initial(8, :) = [1, 0, 0, 0, 1, 1];
x_initial(9, :) = [1, 1, 0, 0, 0, 1];
x_initial(10, :) = [1, 1, 0, 1, 1, 1];

y_experiment = [];
spectral_estimates = [];
M = 2^6 -1 ; 
%now generate all y's
for i = 1:10
    y_experiment(i, :) = generate_prbs(gamma, x_initial(i, :));
    y_period = y_experiment(i, 1:M) ;
    autocorrelation_y_periodic = autocorrelation_periodic(y_period, lags);
    spectral_estimates(i, :) = fft(autocorrelation_y_periodic);
end
averaged_spectral_estimate = mean(spectral_estimates, 1);
T_s = 1;
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
m = size(W_t, 2);

figure(3);
subplot(3, 1, 1);
% plot(0:M-1, autocorrelation_y_periodic);
loglog(W_t, abs(averaged_spectral_estimate(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('spectral of PRBS');
axes.Title.FontSize = 18;

y_experiment = [];
spectral_estimates = [];
% Generate a normaly distributed signal
for i = 1:10
    y_norm = randn(M, 1)*gamma;
    %do the saturation
    for j = 1:M
        if y_norm(j) > gamma
            y_norm(j) = gamma ;
        elseif y_norm(j) < -gamma
            y_norm(j) = -gamma ;
        end
    end

    y_experiment(i, :) = y_norm;
    autocorrelation_y_periodic = autocorrelation_periodic(y_experiment(i, :), lags);
    spectral_estimates(i, :) = fft(autocorrelation_y_periodic);
end
averaged_spectral_estimate = mean(spectral_estimates, 1);
T_s = 1;
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
m = size(W_t, 2);

subplot(3, 1, 2);
loglog(W_t, abs(averaged_spectral_estimate(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('spectral of normal');
axes.Title.FontSize = 18;

% Generate a uniformly distributed signal
y_experiment = [];
spectral_estimates = [];
for i = 1:10
    y_uniform = -gamma + 2*gamma.*rand(M,1);
    y_experiment(i, :) = y_uniform;
    autocorrelation_y_periodic = autocorrelation_periodic(y_experiment(i, :), lags);
    spectral_estimates(i, :) = fft(autocorrelation_y_periodic);
end
averaged_spectral_estimate = mean(spectral_estimates, 1);
T_s = 1;
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
m = size(W_t, 2);

subplot(3, 1, 3);
loglog(W_t, abs(averaged_spectral_estimate(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('spectral of uniform');
axes.Title.FontSize = 18;
