function [p1_R,p1_omega,p1_a,p1_var] = HS2019_SysID_midterm_p1_18936872()
%% Solution for Problem 1
%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data

% Extract Legi from Filename
clear;
name=mfilename;
LegiNumber= name(end-7:end);

p1_U = HS2019_SysID_midterm_p1_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variable p1_U to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Task 1: Calculation of Autocorrelation
%Initialize

clc;
close all;
screensize = get(groot, 'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);

%generate fake data
%generate u(k) periodic

% r = 5;
% M = 1024;
% L = r*M;
% u_period = 2 * randn(1, M);
% u = repmat(u_period, 1, r);
%end fake data

 %calculate period using each record seperatly
 num_records = size(p1_U, 1);
 period_records = zeros(1, num_records);
 for record=1:num_records
    u = p1_U(record, :);
    L = size(u, 2);
    sampling_frequency = 1/0.5;
    t = (0:L - 1)/sampling_frequency;
    %compute autocorrelation
    lags = ((-L/2) + 1):(L/2);
    autocorrelation_u_periodic = autocorrelation_periodic(u, lags);
    [pksh,lcsh] = findpeaks(autocorrelation_u_periodic);
    short = mean(diff(lcsh))/sampling_frequency;
    [pklg,lclg] = findpeaks(autocorrelation_u_periodic, 'MinPeakDistance',ceil(short)*sampling_frequency,'MinPeakheight',70);
    period_records(record) = mean(diff(lclg))/sampling_frequency;
 end
%finally the calculated period is 
period = mean(period_records);
T = period * sampling_frequency;
%plot autocorrelation of the first window
u = p1_U(1, 1:T);
lags = (((-T)/2) + 1):(T/2);
autocorrelation_u_periodic = autocorrelation_periodic(u, lags);

figure(1);
plot(lags, autocorrelation_u_periodic);
axis tight;
axes = gca;
axes.Title.Interpreter = 'latex';
axes.Title.String = 'Input signal autocorrelation of one period';
axes.Title.FontSize = 18;
axes.XAxis.TickLabelInterpreter = 'latex';
axes.XAxis.FontSize = 10;
axes.YAxis.TickLabelInterpreter = 'latex';
axes.YAxis.FontSize = 10;
axes.XLabel.Interpreter = 'latex';
axes.XLabel.String = 'lags $\tau$';
axes.XLabel.FontSize = 14;
axes.YLabel.Interpreter = 'latex';
axes.YLabel.String = '$R_{uu}(T)$';
axes.YLabel.FontSize = 14;

p1_R = autocorrelation_u_periodic;

disp("===========================================================================================================================================");
disp("Problem 1.1");
disp("First we plot the autocorrelation of every record provided and then we identify the peaks which are then used to estimate the period length");
disp("through calculating the distance between the peaks themselves");
fprintf("calculated period is %d\n", period);
disp("===========================================================================================================================================");

%% Task 2: Estimation of signal components



L = 500;
fs = 1;
t = (0:L-1) .* (1/fs);
freq = 100;
u = 37 * cos(((2*pi)/freq)*t);

% find_frequencies_amplitude(u, fs)

u_period = u(1: freq);
U_freq = fft(u_period);
% figure(2);
omega = (0:freq-1) * (2*pi/freq);
idx = find(omega > 0 & omega < pi);
plot(omega(idx), abs(U_freq(idx)), 'linewidth', 2, 'Color', 'red');
% 
% [pksh,lcsh] = findpeaks(U_freq(idx), 'MinPeakheight',0.5*(10^4));
% frequencies = omega(idx);
% frequencies = frequencies(lcsh) .* sampling_frequency;
% amplitude = Pu(idx);
% amplitude = sqrt(amplitude(lcsh));
        

function [frequencies, amplitude] = find_frequencies_amplitude(u, sampling_frequency)
    %calculate the periodogram of w
    N = size(u, 2);
    t = 0:N-1;
    omega = t * (2*pi/N);
    idx = find(omega > 0 & omega < pi);
    Pu = periodogram(u);
    %plotting time
    figure(3);
    plot(omega(idx), Pu(idx), 'linewidth', 2, 'Color', 'red');
    fig = gcf;
    fig.Position = [screenwidth/2, 0, screenwidth/4, screenheight/4];
    axes = gca;
    axes.Title.Interpreter = 'latex';
    [pksh,lcsh] = findpeaks(Pu(idx), 'MinPeakheight',0.5*(10^4));
    frequencies = omega(idx);
    frequencies = frequencies(lcsh) .* sampling_frequency;
    amplitude = Pu(idx);
    amplitude = sqrt((amplitude(lcsh).*4) ./N);
end

samples_per_record = size(p1_U, 2);
number_of_periods_per_record = floor(samples_per_record/T);
remainder = mod(samples_per_record, T);
omega = (0:T-1) * (2*pi/T);
idx = find(omega > 0 & omega < pi);
for record=1:num_records
    u = p1_U(record, :);
    u = reshape(u(remainder+1:samples_per_record), [T, number_of_periods_per_record]);
    %loop around each period calculating the estimates
    for period=1:number_of_periods_per_record
        u_period = u(:, period);
        U_freq = periodogram(u);
        [pksh,lcsh] = findpeaks(U_freq(idx), 'MinPeakheight',0.5*(10^4));
        frequencies = omega(idx);
        frequencies = frequencies(lcsh) .* sampling_frequency;
        amplitude = Pu(idx);
        amplitude = sqrt(amplitude(lcsh));
        
        figure(2);
        plot(omega, U_freq, 'linewidth', 2, 'Color', 'red');
    end
    
    %get fft of each period
 end
% 

fprintf('The frequencies for the fist chunk is\n');
[first_chunk_freq, first_chunk_amps] = find_frequencies_amplitude(p1_U(1, :), sampling_frequency);
disp(first_chunk_freq);
fprintf('The amplitudes for the fist chunk is\n');
disp(first_chunk_amps);
fprintf('The frequencies for the second chunk is\n');
[second_chunk_freq, second_chunk_amps] = find_frequencies_amplitude(p1_U(2, :), sampling_frequency);
disp(second_chunk_freq);
fprintf('The amplitudes for the second chunk is\n');
disp(second_chunk_amps);
fprintf('The frequencies for the third chunk is\n');
[third_chunk_freq, third_chunk_amps] = find_frequencies_amplitude(p1_U(3, :), sampling_frequency);
disp(third_chunk_freq);
fprintf('The amplitudes for the third chunk is\n');
disp(third_chunk_amps);
fprintf('The frequencies for the fourth chunk is\n');
[fourth_chunk_freq, fourth_chunk_amps] = find_frequencies_amplitude(p1_U(4, :), sampling_frequency);
disp(fourth_chunk_freq);
fprintf('The amplitudes for the fourth chunk is\n');
disp(fourth_chunk_amps);


% L = 10000000;
% fs = 8;
% t = (0:L-1) .* (1/fs);
% u = 37 * sin(((2*pi)/4)*t);
% find_frequencies_amplitude(u, fs)


p1_omega = first_chunk_freq;
p1_a = (first_chunk_amps + second_chunk_amps + third_chunk_amps + fourth_chunk_amps)./4;

%% Task 3: Estimation of noise variance
    figure(3);
    N = size(p1_U, 2);
    t = 0:N-1;
    plot(t, p1_U(1, :));
    fig = gcf;
    fig.Position = [0, screenheight, screenwidth/4, screenheight/4];
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Input u signal';

p1_var = zeros(1, 5);
for i=1:5
    %resize
    N = size(p1_U, 2);
    number_of_periods = floor(N/T);
    remainder = mod(N, T);
    reshaped_record = reshape(p1_U(i, remainder+1:N), [number_of_periods, T]);
    %get mean
    figure(4);
    plot(0:T-1, reshaped_record(1, :), 'linewidth', 2, 'Color', 'red');
    hold on;
    plot(0:T-1, reshaped_record(2, :), 'linewidth', 2, 'Color', 'blue');
    hold on;
    plot(0:T-1, reshaped_record(3, :), 'linewidth', 2, 'Color', 'green');
    mean_record = mean(reshaped_record, 1);
    hold on;
    plot(0:T-1, mean_record, 'linewidth', 2, 'Color', 'black');
    fig = gcf;
    fig.Position = [screenwidth/2, 0, screenwidth/4, screenheight/4];
    axes = gca;
    axes.Title.Interpreter = 'latex';
    
    variance_record = mean((reshaped_record - repmat(mean_record, number_of_periods, 1)).^2, 1);
    p1_var(i) = mean(variance_record);
end
p1_var = mean(p1_var, 2);

    function result = autocorrelation_periodic(x, shifts)
        N_x = size(x, 2); %period
        num_shifts = size(shifts, 2);
        x_shifted = zeros(N_x, num_shifts);
        for s=1:num_shifts
            %shift the period signal
            x_shifted(:, s) = circshift(x, shifts(1, s))';
        end
        R_x = 1/N_x * x * x_shifted;
        result = R_x;

    end
    function result = periodogram(e)
            N_periodogram = size(e, 2);
            En = fft(e);
            result = 1/N_periodogram .* (abs(En) .^ 2);
    end
end

