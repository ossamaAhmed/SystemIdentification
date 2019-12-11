clc;
close all;

screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%generate random signal
N = 4096;
e = randn(1, N);

%calculate the periodogram
Pe = periodogram(e);

%plot on log-log scale
w_n = 0:N-1;
w_n = w_n * (2*pi/N);
figure(1);
loglog(w_n, Pe, 'Color', 'blue');
fig = gcf;
fig.Position = [0, 0, screenwidth/4, screenheight/4];

%DT Plant
num = 1;
denom = [1, -0.9, 0.5];
sampling_time = 1;
plant = tf(num, denom, sampling_time);

w = lsim(plant, e);

%calculate the periodogram of w
Pw = periodogram(w');

%calculate the magnitude of the bode plot of the transfer function
[mag,phase,wout] = bode(plant);

% Get magnitude of frequency response of P squared
plant_freq_response = freqresp(plant, w_n);
plant_abs_squared = abs(plant_freq_response(:))' .^2;

%plotting time
figure(2);
loglog(w_n, Pw, 'linewidth', 2, 'Color', 'red');
axis tight;
hold on;
loglog(w_n, plant_abs_squared, 'linewidth', 4, 'Color', 'blue');
fig = gcf;
fig.Position = [screenwidth/2, 0, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';


err = Pw - plant_abs_squared;
figure(3);
semilogx(w_n, err);
axis tight;
fig = gcf;
fig.Position = [0, screenheight/2, screenwidth/4, screenheight/4];
