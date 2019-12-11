% Initialization
clc;
clear all;
close all;
% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%define the plant first
s = tf('s');
Eta_z = 0.1;
W_z = 3;
Eta_p = 0.1;
W_p = 3.5;
Num_G = (s^2 + (2 * Eta_z * W_z * s) + (W_z^2)) * 5000;
Denom_G = (s^2 + (2 * Eta_p * W_p * s) + (W_p^2)) * (s + 50) * (s + 200);
G_s = Num_G / Denom_G;
T_s = 0.02;
G_d = c2d(G_s, T_s, 'zoh');

%deine the controller (PI)
z = tf('z', T_s);
C_d = ((1.25 * z) - 0.75) / (z - 1);
lambda = 0.1;
std_dev = sqrt(lambda);

%------------------------------------------------------
%------------------PART ONE----------------------------
%------------------------------------------------------
%calculate the discrete-time senstivity function
S_d = 1 / (1 + (G_d * C_d));
T_d = 1 - S_d;
cl_tf = feedback(G_d*C_d,-1);
%check if all poles of S_d are inside the unit circle
S_d_poles = pole(S_d);
T_d_poles = pole(T_d);
if abs(S_d_poles) < 1
    disp('checked all poles are inside the unit circle');
else
    disp('Error: not all poles are inside the unit circle');
end

if S_d_poles == T_d_poles
    disp('checked all poles are the same between S_d and T_d');
else
    disp('Error: not all poles are the same between S_d and T_d');
end
%------------------------------------------------------
%------------------PART TWO----------------------------
%------------------------------------------------------

%create prbs signal
max_amplitude = 0.1;
order = 10; % got it from the formula
M = (2 .^ order) - 1;
Band = [0,1];
Range = [-max_amplitude, max_amplitude];
r_period = idinput([M, 1, 1], 'prbs', Band, Range)';
P = 2;
r = repmat(r_period, 1, P);

%get input and output
v = std_dev * randn(1, P*M);
y = lsim(T_d, r)' + lsim(S_d, v)';
u = lsim(S_d, r)' - lsim(S_d * C_d, v)';

%get estimate now
Gest = zeros(1, M);
for i = 2 : P
    u_period = u((i-1)*M+1:i*M);
    y_period = y((i-1)*M+1:i*M);
    U_period = fft(u_period);
    Y_period = fft(y_period);
    Gest = Gest + Y_period ./ U_period;
end
Gest = Gest / (P - 1);

%get real frequency response of the system
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
Gfreq = freqresp(G_d, W_t);
Gfreq = Gfreq(:)';
m = size(W_t, 2);

%magnitude plot
figure(2);
subplot(2,1,1);
loglog(W_t, abs(Gfreq));
hold on;
loglog(W_t, abs(Gest(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Magnitude';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$\left | \hat{G}_{y}(e^{jw}) \right |$';
axes.YLabel.FontSize = 14;

%phase plot
subplot(2,1,2);
semilogx(W_t, phase(Gfreq));
hold on;
semilogx(W_t, phase(Gest(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Phase';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$arg(\hat{G}_{y}(e^{jw}))$ $[rad]$';
axes.YLabel.FontSize = 14;


%------------------------------------------------------
%------------------PART THREE----------------------------
%------------------------------------------------------
%create prbs signal
max_amplitude = 0.1;
order = 10; % got it from the formula
M = (2 .^ order) - 1;
Band = [0,1];
Range = [-max_amplitude, max_amplitude];
r_period = idinput([M, 1, 1], 'prbs', Band, Range)';
R_period = fft(r_period);
P = 10;
r = repmat(r_period, 1, P);

%get input and output
v = std_dev * randn(1, P*M);
y = lsim(T_d, r)' + lsim(S_d, v)';
e = r - y;
S_d_est = zeros(1, M);
for i = 2 : P
    e_period = e((i-1)*M+1:i*M);
    E_period = fft(e_period);
    S_d_est = S_d_est + (E_period ./ R_period);
end
S_d_est = S_d_est / (P - 1);

%get real frequency response of the system
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
Sfreq = freqresp(S_d, W_t);
Sfreq = Sfreq(:)';
m = size(W_t, 2);

%magnitude plot
figure(3);
subplot(2,1,1);
loglog(W_t, abs(Sfreq));
hold on;
loglog(W_t, abs(S_d_est(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Sensitivity Magnitude';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$\left | \hat{G}_{y}(e^{jw}) \right |$';
axes.YLabel.FontSize = 14;

%phase plot
subplot(2,1,2);
semilogx(W_t, phase(Sfreq));
hold on;
semilogx(W_t, phase(S_d_est(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Sensitivity Phase';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$arg(\hat{G}_{y}(e^{jw}))$ $[rad]$';
axes.YLabel.FontSize = 14;

%------------------------------------------------------
%------------------PART Four----------------------------
%------------------------------------------------------
%create prbs signal
max_amplitude = 0.1;
order = 10; % got it from the formula
M = (2 .^ order) - 1;
Band = [0,1];
Range = [-max_amplitude, max_amplitude];
w_period = idinput([M, 1, 1], 'prbs', Band, Range)';
W_period = fft(w_period);
P = 10;
w = repmat(w_period, 1, P);

%get input and output
v = std_dev * randn(1, P*M);
y = lsim(S_d * G_d, w)' + lsim(S_d, v)';
S_d_G_d_est = zeros(1, M);
for i = 2 : P
    y_period = y((i-1)*M+1:i*M);
    Y_period = fft(y_period);
    S_d_G_d_est = S_d_G_d_est + (Y_period ./ W_period);
end
S_d_G_d_est = S_d_G_d_est / (P - 1);

%get real frequency response of the system
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
S_G_freq = freqresp(S_d * G_d, W_t);
S_G_freq = S_G_freq(:)';
m = size(W_t, 2);

%magnitude plot
figure(4);
subplot(2,1,1);
loglog(W_t, abs(S_G_freq));
hold on;
loglog(W_t, abs(S_d_G_d_est(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Plant * Sensitivity Magnitude';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$\left | \hat{G}_{y}(e^{jw}) \right |$';
axes.YLabel.FontSize = 14;

%phase plot
subplot(2,1,2);
semilogx(W_t, phase(S_G_freq));
hold on;
semilogx(W_t, phase(S_d_G_d_est(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Plant * Sensitivity Phase';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$arg(\hat{G}_{y}(e^{jw}))$ $[rad]$';
axes.YLabel.FontSize = 14;

%------------------------------------------------------
%------------------PART Five----------------------------
%------------------------------------------------------
Gest_second = S_d_G_d_est ./ S_d_est;

%magnitude plot
figure(5);
subplot(2,1,1);
loglog(W_t, abs(Gfreq));
hold on;
loglog(W_t, abs(Gest_second(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Magnitude';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$\left | \hat{G}_{y}(e^{jw}) \right |$';
axes.YLabel.FontSize = 14;

%phase plot
subplot(2,1,2);
semilogx(W_t, phase(Gfreq));
hold on;
semilogx(W_t, phase(Gest_second(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Phase';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = '$\omega$ $[\frac{rad}{sec}]$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$arg(\hat{G}_{y}(e^{jw}))$ $[rad]$';
axes.YLabel.FontSize = 14;
