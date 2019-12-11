%Initialize

clc;
close all;
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%generate u(k) periodic

r = 20;
M = 1024;
L = r*M;
u_period = 2 * randn(1, M);
u = repmat(u_period, 1, r);

%generate output data y
G_num = [0.1, 0];
G_denom = [1, -1.7, 0.72];
sampling_time = 1;
G = tf(G_num, G_denom, sampling_time);

H_num = [1.5, -1.5*0.92];
H_denom = [1, -0.5];
sampling_time = 1;
H = tf(H_num, H_denom, sampling_time);

%generate noise 
e = 0.1 * randn(1, L);

y = lsim(G, u)' + lsim(H, e)';

figure(1);
plot(0:L-1, u);
fig = gcf;
fig.Position = [0, screenheight, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Input u signal';

figure(2);
plot(0:L-1, y);
fig = gcf;
fig.Position = [screenwidth/2, screenheight, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Response Y signal';

%compute autocorrelation
lags = ((-L/2) + 1):(L/2);
autocorrelation_u_periodic = autocorrelation_periodic(u_period, lags);
autocorrelation_u_finite = autocorrelation_finite(u, lags);

figure(3);
plot(0:L-1, autocorrelation_u_periodic);
fig = gcf;
fig.Position = [0, screenheight/2, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'periodic autocorrelation u  signal';
axes.XLabel.String = 'lags';

figure(4);
plot(0:L-1, autocorrelation_u_finite);
fig = gcf;
fig.Position = [screenwidth/2, screenheight/2, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'finite autocorrelation u  signal';
axes.XLabel.String = 'lags';

figure(5);
plot(0:2*L-2, xcorr(u));
fig = gcf;
fig.Position = [screenwidth, 0, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'matlab autocorrelation u signal';
axes.XLabel.String = 'lags';

%average system input and output over r-1 periods
u_excluding_transient = u(1, M+1:L);
y_excluding_transient = y(1, M+1:L);
u_avergaed = sum(reshape(u_excluding_transient,[M,r-1]), 2) ./ (r-1);
y_avergaed = sum(reshape(y_excluding_transient,[M,r-1]), 2) ./ (r-1);
u_avergaed = u_avergaed';
y_avergaed = y_avergaed';

%calculate the ETFE now
U_frequency_response = fft(u_avergaed);
Y_frequency_response = fft(y_avergaed);
G_estimated = Y_frequency_response ./ U_frequency_response;
omega = (2*pi/M)*[0:M-1];
idx = find(omega > 0 & omega < pi);

%find the real response
Gfreq = freqresp(G, omega);
Gabs = abs(Gfreq(:))';


figure(6);
loglog(omega(idx), abs(G_estimated(idx)), 'Color', 'blue');
hold on;
loglog(omega(idx), abs(Gabs(idx)), 'Color', 'Red');
fig.Position = [screenwidth, 0, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Estimate VS Real Response';

RMS = rms(Gabs - G_estimated);
fprintf("RMS Error is: %d\n", RMS);


