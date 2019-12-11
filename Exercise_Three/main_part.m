% %Initialize

clc;
close all;
clear;
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%generate u(k) periodic
N = 1024;
u = 10 * randn(1, N);

%generate output data y
G_num = [0.1, 0];
G_denom = conv([1 -1.7 0.72], [1 -0.98 0.9]);
sampling_time = 1;
G = tf(G_num, G_denom, sampling_time);

H_num = [0.5, -0.9*0.5];
H_denom = [1, -0.25];
sampling_time = 1;
H = tf(H_num, H_denom, sampling_time);

%generate noise 
e = randn(1, N);

y = lsim(G, u)' + lsim(H, e)';

figure(1);
plot(0:N-1, u);
fig = gcf;
fig.Position = [0, screenheight, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Input u signal';

figure(2);
plot(0:N-1, y);
fig = gcf;
fig.Position = [screenwidth/2, screenheight, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Response Y signal';

%Now estimate G(Z) via the unsmooted ETFE and compare

U_frequency_response = fft(u);
Y_frequency_response = fft(y);
G_estimated = Y_frequency_response ./ U_frequency_response;
omega = (2*pi/N)*[0:N-1];
idx = find(omega > 0 & omega < pi);

%find the real response
Gfreq = freqresp(G, omega);
Gfreq = Gfreq(:)';



figure(3);
loglog(omega(idx), abs(G_estimated(idx)), 'Color', 'Blue');
hold on;
loglog(omega(idx), abs(Gfreq(idx)), 'Color', 'Green');
% loglog(omega(idx), abs(G_estimated(idx)- Gabs(idx)), 'Color', 'Red');
fig.Position = [screenwidth, 0, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Estimate VS Real Response';

%find the real response of noise tf
Hfreq = freqresp(H, omega);
Hfreq = Hfreq(:)';

figure(4);
loglog(omega(idx), abs(G_estimated(idx) - Gfreq(idx)), 'Color', 'Blue');
hold on;
loglog(omega(idx), abs(Hfreq(idx)), 'Color', 'Red');
fig.Position = [screenwidth, 0, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Real Error';

%split data into 4 chunks
R = 4;
M = N / R;
omega_per_chunk = (2*pi/M)*[0:M-1];
idx_per_chunk = find(omega_per_chunk > 0 & omega_per_chunk < pi);
G_estimated_all = zeros(size(idx_per_chunk, 2), R);
alpha_sum = 0;
for i = 1:R
    %transform 
    U_frequency_response_i = fft(u(1, ((i-1)*M)+1:i*M));
    Y_frequency_response_i = fft(y(1, ((i-1)*M)+1:i*M));
    G_estimated_i = Y_frequency_response_i ./ U_frequency_response_i;
    alpha = periodogram(u(1, ((i-1)*M)+1:i*M));
    alpha_sum = alpha_sum + alpha;
    G_estimated_all(:, i) = (alpha(idx_per_chunk) .* G_estimated_i(idx_per_chunk))';
end
G_smoothed_average = sum(G_estimated_all, 2) ./ alpha_sum(idx_per_chunk)';

figure(5);
loglog(omega(idx), abs(G_estimated(idx)), 'Color', 'Blue');
hold on;
loglog(omega(idx), abs(Gfreq(idx)), 'Color', 'Green');
hold on;
loglog(omega_per_chunk(idx_per_chunk), abs(G_smoothed_average), 'Color', 'Red');
fig.Position = [screenwidth, 0, screenwidth/4, screenheight/4];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'ETFE Estimate VS Real Response';

%plotting the frequency response of the hann window for different gammas
g5 = 5;
g10 = 10;
g50 = 50;
g100 = 100;
[omegag5, WHg5] = WfHann(g5, N);
[omegag10, WHg10] = WfHann(g10, N);
[omegag50, WHg50] = WfHann(g50, N);
[omegag100, WHg100] = WfHann(g100, N);

figure(6);
plot(omegag5, WHg5);
hold on;
plot(omegag10, WHg10);
hold on;
plot(omegag50, WHg50);
hold on;
plot(omegag100, WHg100);
legend('gamma 5','gamma 10', 'gamma 50', 'gamma 100');
fig = gcf;
fig.Position = [0, screenheight/4, screenwidth/2, screenheight/2];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Frequency Domain Hann Windows';

%Now calculate a smoothed ETFE with the window of Hann
[ETFE_Hann_Estimate_gamma_5] = ETFE_Smoothed_Hann(u, y, g100);
figure(7);
loglog(omega(idx), abs(ETFE_Hann_Estimate_gamma_5(idx)));
hold on;
loglog(omega(idx), abs(Gfreq(idx)), 'Color', 'Green');
fig = gcf;
fig.Position = [0, screenheight/4, screenwidth/2, screenheight/2];
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Frequency Domain Hann Windows Estimate';
