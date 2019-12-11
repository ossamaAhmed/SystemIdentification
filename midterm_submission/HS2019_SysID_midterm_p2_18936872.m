function [p2_umin, p2_umax, p2_M, p2_Np, p2_u_etfe, p2_y_etfe, p2_omega, p2_estimate_no_window, p2_estimate_windowed, p2_gamma_best] = HS2019_SysID_midterm_p2_18936872()

    %% Solution for Problem 2
    
    %% General instructions for solution
    % Change the filename of this function, both in the function definition
    % above and in the filename in the folder 

    % Modify your code in the next sections, and return the variables
    % requested.

    % If you skip one part of the problem, return the empty vectors as already
    % provided in the code

    % use the plant-function like this (where p2_u is any input you wish 
    % to apply and p2_y is the corresponding output of the plant):    
    % p2_y = HS2019_SysID_midterm_p2_system_sim(LegiNumber,p2_u);

    % use the validation-function like this (where p2_M is the period length 
    % of p2_u_etfe you found in Task 2):
    % p2_u_cv = HS2019_SysID_midterm_p2_validation(p2_M);
    % Extract Legi from Filename
    clear;
    name = mfilename;
    LegiNumber = str2double(name(end-7:end));
    
    %% Task 1: Estimate plant input saturation 
    clc;
    close all;
    screensize = get(groot, 'ScreenSize');
    screenwidth = screensize(3);
    screenheight = screensize(4);

    %create an input signal formed of step functions
    u_range = -5:0.1:5;
    T = 90;
    u = repmat(u_range, T, 1)';
    u = reshape(u.',1,[]);
    %get corresponding output
    y = HS2019_SysID_midterm_p2_system_sim(LegiNumber,u);
    %smooth y output by averaging per unit step
    y_reshaped = reshape(y, T, size(u_range, 2));
    y_smoothed = mean(y_reshaped, 1);
    t = 0:(T*size(u_range, 2))-1;
    t_smoothed = 0:T:(T*size(u_range, 2))-1;
   
    y_gradient = gradient(y_smoothed);
    indicies_for_neg_gradient_lower = find(y_gradient(1:ceil(size(y_gradient, 2)/2))<0);
    indicies_for_neg_gradient_upper = find(y_gradient(ceil(size(y_gradient, 2)/2):size(y_gradient, 2))<0);
    
    p2_umin = u_range(indicies_for_neg_gradient_lower(end) + 2); %rounded to first decimal
    p2_umax = u_range(indicies_for_neg_gradient_upper(1) + ceil(size(y_gradient, 2)/2) - 4); %rounded to first decimal
    disp("===========================================================================================================================================");
    disp("Problem 2.1");
    disp("The saturation limits were estimated by generating a unit step signal corresponding to the specified true limits [-5, 5]");
    disp("where the same input will be asserted for 90s, afterwards this step signal (stairs) will be be passed through the system to get a response");
    disp("After getting the response, we get the mean response per unit value in the range and then get the gradient of this mean to see where is");
    disp("the saturation level by checking where the signal plateau's and becomes noisy around a certain level");
    fprintf("u minimum is %.1f and u maximum is %.1f\n", p2_umin, p2_umax);
    disp("===========================================================================================================================================");
    
    %% Task 2: Design correct input
    order = 9; % got it from the formula
    M = (2 .^ order) - 1;
    Band = [0,1];
    %decide on which input range to have it symmetric
    u_limit = abs(min(p2_umin, p2_umax));
    Range = [-u_limit, u_limit];
    u_period = idinput([M, 1, 1], 'prbs', Band, Range)';
    %generate one period to calculate the expected variance of the estimate
    expected_estimate_variance = mean(0.25 ./ periodogram(u_period));
    r = ceil(expected_estimate_variance/0.04) + 1;
    L = r*M;
    u = repmat(u_period, 1, r);

    y = HS2019_SysID_midterm_p2_system_sim(LegiNumber,u);
    p2_M       = M; % integer scalar, length of one period in p2_u_etfe
    p2_Np      = r; % integer scalar, number of periods in p2_u_etfe
    p2_u_etfe  = u; % vector, your designed input signal used for the ETFE
    p2_y_etfe  = y'; % vector, the output from the plant, used for the ETFE
    disp("===========================================================================================================================================");
    disp("Problem 2.2");
    disp("To first estimate the period length, we know that the period length KT determines the frequency resolution 2pi/T, since we know the resolution");
    disp("we are trying to reach which is pi/200 then the period length has to be bigger than 400, therefore we choose an order of 9 giving us a period");
    disp("of length 511, for the number of periods, we know that the variance of our estimate(1 period) is spectogram_of_noise/periodgram_of_input, and since");
    disp("our noise is a white noise with zero mean then the spectogram of noise equals to its variance 0.25, once we get the variance of our estimate ");
    disp("we know that variance of our estimate(r periods) = variance_of_estimate_one_period/R, which we leverage to get the number of periods needed");
    disp("===========================================================================================================================================");

    %% Task 3: Compute unsmoothed ETFE
    omega_per_chunk = (2*pi/p2_M)*[0:p2_M-1];
    idx_per_chunk = find(omega_per_chunk > 0 & omega_per_chunk < pi);
    G_estimated_all = zeros(size(idx_per_chunk, 2), p2_Np-1);
    alpha_sum = 0;
    for i = 2:p2_Np
        %transform 
        U_frequency_response_i = fft(p2_u_etfe(1, ((i-1)*p2_M)+1:i*p2_M));
        Y_frequency_response_i = fft(p2_y_etfe(1, ((i-1)*p2_M)+1:i*p2_M));
        G_estimated_i = Y_frequency_response_i ./ U_frequency_response_i;
        alpha = periodogram(p2_u_etfe(1, ((i-1)*p2_M)+1:i*p2_M));
        alpha_sum = alpha_sum + alpha;
        G_estimated_all(:, i-1) = (alpha(idx_per_chunk) .* G_estimated_i(idx_per_chunk))';
    end
    G_un_smoothed_average = sum(G_estimated_all, 2) ./ alpha_sum(idx_per_chunk)';
    figure(1)
    % Magnitude plot
    subplot(2,1,1);
    loglog(omega_per_chunk(idx_per_chunk), abs(G_un_smoothed_average), '-o');
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Unsmoothed ETFE Magnitude';
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
    % Phase plot
    subplot(2,1,2);
    semilogx(omega_per_chunk(idx_per_chunk), angle(G_un_smoothed_average), '-o');
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Unsmoothed ETFE Phase';
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
    p2_omega               = omega_per_chunk(idx_per_chunk); % vector, equally spaced frequecies in (0,pi)(not inclusive) 
    p2_estimate_no_window  = G_un_smoothed_average; % vector, ETFE estimate unsmoothed
    
    disp("===========================================================================================================================================");
    disp("Problem 2.3");
    disp("We follow the same way to estimate the unsmoothed ETFE as in the lectures,  we do this by division of the fft(y)/fft(u) and we weight");
    disp("the averaging process by the inverse of the variance, since we know the signal is periodic we do the averaging on each period estimate,")
    disp("===========================================================================================================================================");
    
    %% Task 4: Compute smoothed ETFE with best Hann window width
    %perform cross validation
    omegas = (2*pi/p2_M)*[0:p2_M-1];
    picked_indicies = find(omegas > 0 & omegas < pi);
    p2_u_cv = HS2019_SysID_midterm_p2_validation(p2_M);
    actual_output = HS2019_SysID_midterm_p2_system_sim(LegiNumber, p2_u_cv);
    p2_u_cv = p2_u_cv';
    actual_output = actual_output';
    %get the average frequency response of the signal
    p2_u_cv_freq = fft(p2_u_cv);
    p2_u_cv_freq = reshape(p2_u_cv_freq, p2_M, 6);
    p2_u_cv_freq = mean(p2_u_cv_freq, 2);
    actual_output_freq = fft(actual_output);
    actual_output_freq = reshape(actual_output_freq, p2_M, 6);
    actual_output_freq = mean(actual_output_freq, 2);
    gammas = 50:50:500;
    num_gammas = size(gammas, 2);
    gamma_errors = zeros(num_gammas, 1);
    for g=1:num_gammas
        current_gamma = gammas(g);
        G_estimated_all_smoothed = zeros(p2_M, p2_Np-1);
        alpha_sum = 0;
        for i = 2:p2_Np
            [G_estimated_i] = ETFE_Smoothed_Hann(p2_u_etfe(1, ((i-1)*p2_M)+1:i*p2_M), p2_y_etfe(1, ((i-1)*p2_M)+1:i*p2_M), current_gamma);
            alpha = periodogram(p2_u_etfe(1, ((i-1)*p2_M)+1:i*p2_M));
            alpha_sum = alpha_sum + alpha;
            G_estimated_all_smoothed(:, i-1) = (alpha .* G_estimated_i)';
        end
        G_smoothed_average = sum(G_estimated_all_smoothed, 2) ./ alpha_sum';
        estimated_output_freq = G_smoothed_average .* p2_u_cv_freq;
        error = actual_output_freq(picked_indicies(1:40)) - estimated_output_freq(picked_indicies(1:40));
        gamma_errors(g, :) = sqrt(sum(abs(error).^2));
    end
    figure(2);
    plot(gammas, gamma_errors, '-o');
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Error Equation';
    axes.Title.FontSize = 18;
    axes.XAxis.TickLabelInterpreter = 'latex';
    axes.XAxis.FontSize = 10;
    axes.YAxis.TickLabelInterpreter = 'latex';
    axes.YAxis.FontSize = 10;
    axes.XLabel.Interpreter = 'latex';
    axes.XLabel.String = '$\gamma$';
    axes.XLabel.FontSize = 14;
    axes.YLabel.Interpreter = 'latex';
    axes.YLabel.String = '2-norm error between true Y and estimated Y';
    axes.YLabel.FontSize = 14;
    
    min_error = min(gamma_errors);
    [best_gamma_idx] = find(gamma_errors==min_error);
    best_gamma = gammas(best_gamma_idx);
    G_estimated_all_smoothed = zeros(p2_M, p2_Np-1);
    alpha_sum = 0;
    for i = 2:p2_Np
        [G_estimated_i] = ETFE_Smoothed_Hann(p2_u_etfe(1, ((i-1)*p2_M)+1:i*p2_M), p2_y_etfe(1, ((i-1)*p2_M)+1:i*p2_M), best_gamma);
        alpha = periodogram(p2_u_etfe(1, ((i-1)*p2_M)+1:i*p2_M));
        alpha_sum = alpha_sum + alpha;
        G_estimated_all_smoothed(:, i-1) = (alpha .* G_estimated_i)';
    end
    G_smoothed_average = sum(G_estimated_all_smoothed, 2) ./ alpha_sum';
    
    p2_estimate_windowed = G_smoothed_average(picked_indicies); %vector, Hann window smoothed estimate using window width p2_gamma_best
    p2_gamma_best        = best_gamma; %scalar integer, best Hann window width
    disp("===========================================================================================================================================");
    disp("Problem 2.4");
    disp("To calculate the error, I only considered frequencies around the resonant frequency");
    disp("1- The value of the error equation changes with varying gamma since windowing adds a bias to reduce the variance of the estimate, which in");
    disp("turn reduces the error till a certain point (bias-variance trade off)");
    disp("2- If we take more periods then the variance of our estimate will decrease and probably we would need a smaller window/ gamma to reach the mimimum error")
    disp("===========================================================================================================================================");
    function result = periodogram(e)
        N = size(e, 2);
        En = fft(e);
        result = 1/N .* (abs(En) .^ 2);
    end

    function [G_result] = ETFE_Smoothed_Hann(u, y, gamma)
        U_freq = fft(u);
        Y_freq = fft(y);
        G_estimated = Y_freq ./ U_freq;
        N = size(G_estimated, 2);
        [omega, Wg] = WfHann(gamma, N);
        zidx = find(omega==0);
        omega = [omega(zidx:N); omega(1:zidx-1)];
        Wg = [Wg(zidx:N), Wg(1:zidx-1)];
        a = U_freq.*conj(U_freq);
        G_result = zeros(1, N);
        for wn =1:N
            Wnorm = 0;
            for xi = 1:N
                widx = mod(xi-wn, N)+1;
                G_result(wn) = G_result(wn) + (Wg(widx) .* G_estimated(xi) .* a(xi));
                Wnorm = Wnorm + (Wg(widx) .* a(xi));
            end
            G_result(wn) = G_result(wn) ./ Wnorm;
        end
    end
    
end








