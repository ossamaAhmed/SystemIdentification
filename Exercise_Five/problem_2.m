% Initialization
clc;
clear all;
close all;

% Screen size used to place plots
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

U = [];
U_o = [];
U_d = [];
for experiment=1:10
    M = 2^8;
    betta = 0.1;
    alpha_d = 1/(2^9);
    u = randn(M, 1);
    u_o = u + betta;
    u_d = u;
    for k =1:M
        u_d(k) = u_d(k) + (alpha_d*(k - 2^7));
    end
    U(experiment, :) = fft(u);
    U_o(experiment, :) = fft(u_o);
    U_d(experiment, :) = fft(u_d);
end
averaged_U = mean(U, 1);
averaged_U_o = mean(U_o, 1);
averaged_U_d = mean(U_d, 1);
T_s = 1;
W = 0 : 2*pi/M : pi;
W_t = W / T_s;
m = size(W_t, 2);
subplot(3, 1, 1);
loglog(W_t, abs(averaged_U(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('spectral of original');
axes.Title.FontSize = 18;
subplot(3, 1, 2);
loglog(W_t, abs(averaged_U_o(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('spectral of offset');
axes.Title.FontSize = 18;
subplot(3, 1, 3);
loglog(W_t, abs(averaged_U_d(1:m)));
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = strcat('spectral of drift');
axes.Title.FontSize = 18;

