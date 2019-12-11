L = 600;
fs = 1;
t = (0:L-1) .* (1/fs);
freq = 150;
u = 1289 .* cos(((2*pi)/freq)*t);

% find_frequencies_amplitude(u, fs)

u_period = u(1: freq);
disp(max(u_period));
U_freq = fft(u_period);
% figure(2);
omega_all = (0:L-1) * (2*pi/L);
idx_omega_all = find(omega_all > 0 & omega_all < pi);
omega = (0:freq-1) * (2*pi/freq);
idx = find(omega > 0 & omega < pi);
plot(omega(idx), abs(U_freq(idx)), 'linewidth', 2, 'Color', 'red');
% plot(omega_all(idx_omega_all), abs(U_freq(idx_omega_all)), 'linewidth', 2, 'Color', 'red');