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
clc;
close all;
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
 sampling_frequency = 1/0.5;
 for record=1:num_records
    u = p1_U(record, :);
    L = size(u, 2);
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

samples_per_record = size(p1_U, 2);
number_of_periods_per_record = floor(samples_per_record/T);
remainder = mod(samples_per_record, T);
omega = (0:T-1) * (2*pi/T);
idx = find(omega > 0 & omega < pi);
fft_periods = [];
psd_periods = [];
for record=1:num_records
    u_i = p1_U(record, :);
    u_i = reshape(u_i(remainder+1:samples_per_record), [T, number_of_periods_per_record]);
    frequencies_periods = zeros(number_of_periods_per_record-1, 4);
    amplitudes_periods = zeros(number_of_periods_per_record-1, 4);
    for p=2:number_of_periods_per_record
        U_freq_i = fft(u_i(:, p));
        fft_periods((record - 1)*(number_of_periods_per_record - 1) + (p-1), :) = U_freq_i;
        Psd_U((record - 1)*(number_of_periods_per_record - 1) + (p-1), :) = 1/T * abs(U_freq_i).^2;
        %get only frequencies I care about
        U_freq_i = U_freq_i(idx);
        [pksh,lcsh] = findpeaks(abs(U_freq_i), 'MinPeakheight', 320);
        frequencies = omega(idx);
        frequencies_periods(p-1, :) = frequencies(lcsh) .* sampling_frequency;
        amplitudes_periods(p-1, :) = 2 .* (abs(U_freq_i(lcsh)) ./T);
    end
    frequencies_records(record, :) = mean(frequencies_periods, 1);
    amplitudes_records(record, :) = mean(amplitudes_periods, 1);
end
figure(11); 
Ts = 1/sampling_frequency;
omegas = 2*pi/(number_of_periods_per_record*T)/Ts*[0:number_of_periods_per_record*T-1];

% Compute average spectrum
fft_periods = [];
psd_periods = [];
for record = 1:num_records
    fft_periods(record,:) = fft(p1_U(record,1:number_of_periods_per_record*T));
    psd_periods(record,:) = 1/(number_of_periods_per_record*T) * abs(fft_periods(record,1:number_of_periods_per_record*T)).^2;
end
Psd_u = mean(psd_periods,1);

% Plot average spectrum to identify peaks
figure(20);
plot(omegas,Psd_u, '*')
xlabel('\omega [rad/sec]')
title('Average PSD')

idx_greater = find(Psd_u > 100 );

%average across records too
final_frequencies = mean(frequencies_records, 1);
final_amplitudes = mean(amplitudes_records, 1);
%sort them now according to amplitude
[final_amplitudes, sorting_idx] = sort(final_amplitudes);
final_frequencies = final_frequencies(sorting_idx); %convert to radians/sec
p1_omega = final_frequencies;
p1_a = final_amplitudes;

disp("===========================================================================================================================================");
disp("Problem 1.2");
disp("First we get the fft of each of the periods provided in all the records while discarding the remainder at the beginning and the first period in each");
disp("record for transient effect, after that we get the peaks and the corresponding frequencies as well as the magnitude which is then used to calculate the amplitude");
disp("which basically is obtained by dividing the magnitude by the period length multiplied by 2 since its a cosine wave, p1_a=2*(peaks(abs(fft(u)))/T)");
disp("===========================================================================================================================================");

%% Task 3: Estimation of noise variance

    %construct a new signal one period
    test_t = (0:T-1) .* (1/sampling_frequency);
    signal_components = zeros(size(p1_omega, 2), size(test_t, 2));
    for i=1:size(p1_omega, 2)
        signal_components(i, :) = p1_a(i) .* cos(p1_omega(i).*test_t);
    end
    u_hat = sum(signal_components, 1);
    
    
    shifts = ((-T/2) + 1):(T/2);
    records_variance = size(1, num_records);
    records_mean = size(1, num_records);
    for record=1:num_records
        u_i = p1_U(record, :);
        u_i = reshape(u_i(remainder+1:samples_per_record), [T, number_of_periods_per_record]);
        %get auto correlation between two signals 
        cross_correlation = crosscorrelation_finite(u_i(:, 5)', u_hat, shifts);
%         figure(4);
%         plot(shifts, cross_correlation);
        max_cross = max(cross_correlation);
        max_cross_idx = find(cross_correlation == max_cross);
        shift_lag = shifts(max_cross_idx);
        shift_u_hat = circshift(u_hat, shift_lag);
        shift_u_hat = repmat(shift_u_hat, number_of_periods_per_record, 1)';
        error = shift_u_hat - u_i;
        records_variance(record) = mean(var(error, 1));
        records_mean(record) = mean(mean(error, 1));
    end
p1_var = mean(records_variance, 2);
disp("===========================================================================================================================================");
disp("Problem 1.3");
disp("I tried two different methods and both yielded same solution. 1- reconstruct the signal using the estimated quantities and then shift the signal based");
disp("on the cross correlation between the constructed one period signal and one period of the measured signal. once we align the two signals, we calculate the errors");
disp("and from there we calculate the variance of the error, we do this for every period availabe. 2- we measure the varaince across periods by caclulating a mean ");
disp("period and then subtracting the periods from this mean period to calculate the variance afterwards. Variance between periods should come from the noise only.");
disp("===========================================================================================================================================");

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


function result = crosscorrelation_finite(y, x, shifts)
    N = size(x, 2); %period
    num_shifts = size(shifts, 2);
    x_shifted = zeros(N, num_shifts);
    for j=1:num_shifts
        %shift the period signal
        s = shifts(1,j);
        if s >= 0
            x_shifted(s+1:N,j) = x(1,1:N-s)';
        else
            x_shifted(1:N+s,j) = x(1,-s+1:N);
        end
    end
    R_y_x = y * x_shifted;
    result = R_y_x;
end
end

