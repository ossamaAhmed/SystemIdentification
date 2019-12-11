function [p1_R,p1_omega,p1_a,p1_var] = HS2019_SysID_midterm_p1_15819279()
%% Solution for Problem 1
%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data
clc
% Extract Legi from Filename
% name = mfilename;
% LegiNumber = name(end-7:end);
LegiNumber = '15819289';
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
clc
close all
u = p1_U';
N = size(u,1);
% U = fft(u);
% filter low power frequencies
% freq_threshold = 1000;
% filter = U > freq_threshold;
% U_filtered = U.*filter;
% u_filtered = ifft(U_filtered);
% [R,lags] = xcorr(u);

% Estimate period from autocorrelation of first signal
R = autocorrelation_finite_energy(u(:,1));
plot(R)
T = 2073 - 1512 + 1; % by visual inspection
u_T = u(1:T,:);
p1_R = autocorrelation_periodic(u_T);
plot(p1_R);

%% Task 2: Estimation of signal components
clc
close all

% Estimation of omega from the fft highest amplitude signals.
U_T = abs(fft(u_T));
history = zeros(1,N);
for M = 1:N
    [~,~,u_a] = extract_coefficients_1(u_T,U_T,T,M);
    history(M) = norm(u_a-u_T(:,1));
end
figure
plot(history)
M_opt = 562; % visual inspection

%% plot the loss function to find the optimal M.

[p1_a,p1_omega,u_a] = extract_coefficients_1(u_T,U_T,T,M_opt);
figure
plot(u_a(1:50))
hold on
plot(u_T(1:50,1))
legend('artificial','real')


%% Task 3: Estimation of noise variance
disp('EXPLAIN PROCEDURE')
difference = u_T(:,1)-u_a;
p1_var = var(difference);
disp(p1_var)

end


function R = autocorrelation_finite_energy(x)
N = length(x);
R = zeros(1,N);
for tau = 0:N-1
    s = zeros(1,N);
    for k = tau:N-1
        s(k-tau+1) = x(k+1)*x(k-tau+1);
    end
    R(tau+1) = sum(s);
end
R = [flip(R), R];
end


function [R,lag] = autocorrelation_periodic(u)
T = length(u);
lag = -T/2+1:T/2;
R = zeros(1,T);
for i=1:T
    if lag(i) < 0
        for j = 1-lag(i):T
            R(i) = R(i) + u(j)*u(j+lag(i))/(T-abs(lag(i)));
        end
    else
        for j = 1:T-lag(i)
            R(i) = R(i) + u(j)*u(j+lag(i))/(T-abs(lag(i)));
        end
    end
    
end
end

function [amps, freqs] = find_M_peaks(U,M)
% add sampling frequency to output freqs in Hz.
n = size(U,2);
amps = zeros(n,M);
freqs = zeros(n,M);
for i=1:n
    for j=1:M
        [amp,freq] = max(U(:,i));
        U(freq,i) = 0;
        amps(i,j) = amp;
        freqs(i,j)= freq;
    end
end
end


function [p1_a,p1_omega,u_a] = extract_coefficients_1(u_T,U_T,T,M)

[~,freqs] = find_M_peaks(U_T,M);
Ts = 0.5;
freqs = freqs/(T*Ts); % convert to Hertz
p1_omega = zeros(1,M);
for i=1:M
    p1_omega(i) = 2*pi*mean(freqs(:,i));
end

% Least Squares.
A = zeros(T,M);
for i = 1:T
    for j = 1:M
        A(i,j) = cos(p1_omega(j)*i);
    end
end
p1_a = pinv(A)*u_T(:,1);

% Compare signals.
u_a = zeros(T,1);
for t = 1:T
    s = zeros(1,M);
    for k = 1:M
        s(k) = p1_a(k)*cos(p1_omega(k)*t);
    end
    u_a(t) = sum(s);
end

end



