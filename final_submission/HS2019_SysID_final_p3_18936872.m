function [p3_b_ML,p3_b_MAP,p3_cv_error,p3_prior_best] = HS2019_SysID_final_p3_18936872()

    %% Solution for Problem 3
    screensize = get(groot, 'ScreenSize');
    screenwidth = screensize(3);
    screenheight = screensize(4);
    figure_num = 1;

    %% General instructions for solution
    % Change the filename of this function, both in the function definition
    % above and in the filename in the folder

    % Modify your code in the next sections, and return the variables
    % requested.

    % If you skip one part of the problem, return the null variables as already
    % provided in the code

    % Extract Legi from Filename
    name = mfilename;
    LegiNumber = str2double(name(end-7:end));
    
    % Obtain experiment data
    [p3_u,p3_y,p3_u_cv,p3_y_cv] = HS2019_SysID_final_p3_GenerateData(LegiNumber);

    %% Task 1: Maximum likelihood estimate
    %-----------------------------------------------------------------------
    %inspecting the output and input first
    %-----------------------------------------------------------------------
    fprintf('------------------------------------------------------------------\n');
    fprintf('--------------------TASK ONE--------------------------------------\n');
    fprintf('------------------------------------------------------------------\n');
    fprintf('The regressor has the form y(k|theta) = [u(k-1), u(k-2), u(k-3)...,u(k-8)][b_1;b_2;b_3...;b_8]\n');
    fprintf('Y(theta) = Phi * theta, the ML estimate maximizes the log likelihood of observaing the data which meand minimizing -log(P(Y|theta))\n');
    fprintf('This is also equivilant to minimizing -log(P(e=Y-Phi*theta|theta)) the probbaility of observing the error vector\n');
    fprintf('so P(Y|theta) = ((1/sqrt(0.5*PI))^N*exp(-1/0.5 * sum_k=1_to_N(e(k, theta)^2)))\n');
    fprintf('if we take the derivative and set it to zero, we get 4*sum_k=1_N((Y(K)- PHI*THETA)*PHI) = 0\n');
    fprintf('Which after rearranging we end up getting the least square estimate theta = inv((Phi^T)Phi)(Phi^T)Y\n');
    fprintf('Which was expected since the ML estimate for isotropic covariance, is equivilant to the least square estimate\n');
%     figure(figure_num);
%     figure_num = figure_num + 1;
%     fig = gcf;
%     fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
%     plot(p3_u);
%     hold on;
%     plot(p3_y);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'Input vs. Output';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     lgd = legend('$u(k)$', '$y(k)$');
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = 12;
%     lgd.Location = 'southwest';
    %-----------------------------------------------------------------------
    %Forming the Phi Matrix
    %-----------------------------------------------------------------------
    N = length(p3_u);
    num_params = 8;
    Phi = zeros(N, num_params);
    Phi(2,:) = [p3_u(1) 0 0 0 0 0 0 0];
    Phi(3,:) = [p3_u(2) p3_u(1) 0 0 0 0 0 0];
    Phi(4,:) = [p3_u(3) p3_u(2) p3_u(1) 0 0 0 0 0];
    Phi(5,:) = [p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0 0 0 0];
    Phi(6,:) = [p3_u(5) p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0 0 0];
    Phi(7,:) = [p3_u(6) p3_u(5) p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0 0];
    Phi(8,:) = [p3_u(7) p3_u(6) p3_u(5) p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0];
    for k = num_params+1 : N
        Phi(k,:) = [p3_u(k-1) p3_u(k-2) p3_u(k-3) p3_u(k-4) p3_u(k-5) p3_u(k-6) p3_u(k-7) p3_u(k-8)];
    end
    p2_theta_ML = (Phi.'*Phi) \ (Phi.'*p3_y);
    
    %-----------------------------------------------------------------------
    %Self Validation
    %-----------------------------------------------------------------------
    %cutting off the first 8 elements
    y_estimated = Phi * p2_theta_ML;
    estimation_error = p3_y - y_estimated;
%     fprintf('Figure %d shows model output with estimated thetaML, as well as model output error.\n\n', figure_num);
%     figure(figure_num);
%     figure_num = figure_num + 1;
%     fig = gcf;
%     fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
%     subplot(2,1,1);
%     plot(p3_y);
%     hold on;
%     plot(y_estimated);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'Maximum Likelihood Actual vs. Model Output';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     lgd = legend('$y(k)$', '$y_{\hat{\theta}_{ML}}(k)$');
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = 12;
%     lgd.Location = 'northwest';
%     subplot(2,1,2);
%     plot(estimation_error);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'Maximum Likelihood Model Output Error';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     axes.YLabel.Interpreter = 'latex';
%     axes.YLabel.String = '$y(k) - y_{\hat{\theta}_{ML}}(k)$';
%     axes.YLabel.FontSize = 14;

%     estimation_error_var = var(estimation_error(:), 1);
%     fprintf('Estimated Model ouput error variance is: %.4f.\n', estimation_error_var);
%     fprintf('Estimated Model ouput error is: %.4f.\n', mean(estimation_error(:).^2));
    
    %-----------------------------------------------------------------------
    %Cross Validation
    %-----------------------------------------------------------------------
    Phi = zeros(N, num_params);
    Phi(2,:) = [p3_u_cv(1) 0 0 0 0 0 0 0];
    Phi(3,:) = [p3_u_cv(2) p3_u_cv(1) 0 0 0 0 0 0];
    Phi(4,:) = [p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0 0 0 0];
    Phi(5,:) = [p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0 0 0];
    Phi(6,:) = [p3_u_cv(5) p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0 0];
    Phi(7,:) = [p3_u_cv(6) p3_u_cv(5) p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0];
    Phi(8,:) = [p3_u_cv(7) p3_u_cv(6) p3_u_cv(5) p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0];
    for k = num_params+1 : N
        Phi(k,:) = [p3_u_cv(k-1) p3_u_cv(k-2) p3_u_cv(k-3) p3_u_cv(k-4) p3_u_cv(k-5) p3_u_cv(k-6) p3_u_cv(k-7) p3_u_cv(k-8)];
    end
    y_estimated = Phi * p2_theta_ML;
    estimation_error = p3_y_cv - y_estimated;
%     fprintf('Figure %d shows model cross validation output with estimated thetaML, as well as model cross validation output error.\n\n', figure_num);
%     figure(figure_num);
%     figure_num = figure_num + 1;
%     fig = gcf;
%     fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
%     subplot(2,1,1);
%     plot(p3_y_cv);
%     hold on;
%     plot(y_estimated);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'Maximum Likelihood Actual vs. Model Output cross validation';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     lgd = legend('$y(k)$', '$y_{\hat{\theta}_{ML}}(k)$');
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = 12;
%     lgd.Location = 'northwest';
%     subplot(2,1,2);
%     plot(estimation_error);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'Maximum Likelihood Model Output Error';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     axes.YLabel.Interpreter = 'latex';
%     axes.YLabel.String = '$y(k) - y_{\hat{\theta}_{ML}}(k)$';
%     axes.YLabel.FontSize = 14;
% 
%     estimation_error_var = var(estimation_error(:), 1);
%     fprintf('cross validation Estimated Model ouput error variance is: %.4f.\n', estimation_error_var);
%     fprintf('cross validation Estimated Model ouput error is: %.4f.\n', mean(estimation_error(:).^2));

    p3_b_ML = p2_theta_ML;   % vector of dimension 8x1


    %% Task 2: Maximum a posteriori estimates
    fprintf('------------------------------------------------------------------\n');
    fprintf('--------------------TASK TWO--------------------------------------\n');
    fprintf('------------------------------------------------------------------\n');
    fprintf('Since p(theta) = (1/sqrt((2*pi)^n * det(S)))*exo(-1/2 *(theta^T)*(S^-1)*theta)\n');
    fprintf('Since p(y | theta) = ((1/sqrt(0.5*PI))^N*exp(-1/0.5 * sum_k=1_to_N(e(k, theta)^2))) as stated in the last task\n');
    fprintf('Therefore using bayes theorm we get the posterior proportional to p(theta)*p(y|theta)\n');
    fprintf('we get the following equation to minimize argmin_theta -1/0.5 * sum_k=1_to_N((y(k) - Phi*theta)^2) - 1/2 * theta^T * S^-1 * theta\n');
    fprintf('we get the following after reorganizing argmin_theta 1/0.25 (2norm(y-Phi*theta))^2 + theta^T * S^-1 * theta\n');
    fprintf('by taking the derivative and setting it to zero we get the following equation theta = (Phi^T * Phi + (1/4)* S^-1)^-1 * (Phi^T * Y)\n');
    
    fprintf('The MAP estimate is biased towards the prior distribution given, if its close to the real distribution than the MAP will give a better estimate than the straighforward ML estimate\n');
    fprintf('However if we get a bad prior then it will throw our estimate off.\n');
    fprintf('For using S1, S2 and S3 priors, it assumes that there is no  correlation between the coefficients and that they are independent of each other.\n');
    fprintf('For S2 and S3, the priors have a decreasing variance corresponding to the later inputs which means that our model should care less about the very old inputs for instance u(k-8)\n');
    fprintf('Therefore for these priors, the older the input the less variance it has around the zero mean.\n');
    fprintf('For S4 and S5, the covariance matrices has non diagonal terms which means that the choice of any parameter will affect the choice of the rest of the parameters.\n');
    Phi = zeros(N, num_params);
    Phi(2,:) = [p3_u(1) 0 0 0 0 0 0 0];
    Phi(3,:) = [p3_u(2) p3_u(1) 0 0 0 0 0 0];
    Phi(4,:) = [p3_u(3) p3_u(2) p3_u(1) 0 0 0 0 0];
    Phi(5,:) = [p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0 0 0 0];
    Phi(6,:) = [p3_u(5) p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0 0 0];
    Phi(7,:) = [p3_u(6) p3_u(5) p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0 0];
    Phi(8,:) = [p3_u(7) p3_u(6) p3_u(5) p3_u(4) p3_u(3) p3_u(2) p3_u(1) 0];
    for k = num_params+1 : N
        Phi(k,:) = [p3_u(k-1) p3_u(k-2) p3_u(k-3) p3_u(k-4) p3_u(k-5) p3_u(k-6) p3_u(k-7) p3_u(k-8)];
    end
    %form S covariance matrix for the identity matrix
    S_identity = eye(num_params);
    theta_estimated_s_1 = (Phi'*Phi + ((1/4).*inv(S_identity))) \ (Phi'*p3_y);
    
    %form S covariance matrix for the diagonal matrix 1
    S_diagonal_1 = zeros(num_params);
    for i=1:num_params
        for j=1:num_params
            if i==j
                S_diagonal_1(i, j) = 0.8^i;
            end
        end
    end
    theta_estimated_s_2 = (Phi'*Phi + ((1/4).*inv(S_diagonal_1))) \ (Phi'*p3_y);
    
    %form S covariance matrix for the diagonal matrix 2
    S_diagonal_2 = zeros(num_params);
    for i=1:num_params
        for j=1:num_params
            if i==j
                S_diagonal_2(i, j) = 0.5^i;
            end
        end
    end
    theta_estimated_s_3 = (Phi'*Phi + ((1/4).*inv(S_diagonal_2))) \ (Phi'*p3_y);
    
     %form S covariance matrix for the full matrix 1
    S_full_matrix_1 = zeros(num_params);
    for i=1:num_params
        for j=1:num_params
            S_full_matrix_1(i, j) = 0.8^max(i, j);
        end
    end
    theta_estimated_s_4 = (Phi'*Phi + ((1/4).*inv(S_full_matrix_1))) \ (Phi'*p3_y);
    
    %form S covariance matrix for the full matrix 2
    S_full_matrix_2 = zeros(num_params);
    for i=1:num_params
        for j=1:num_params
            S_full_matrix_2(i, j) = 0.5^max(i, j);
        end
    end
    theta_estimated_s_5 = (Phi'*Phi + ((1/4).*inv(S_full_matrix_2))) \ (Phi'*p3_y);
    p3_b_MAP = [theta_estimated_s_1, theta_estimated_s_2, theta_estimated_s_3, theta_estimated_s_4, theta_estimated_s_5];
    %-----------------------------------------------------------------------
    %Self Validation
    %-----------------------------------------------------------------------
    y_estimated = Phi * theta_estimated_s_1;
    estimation_error = p3_y - y_estimated;
%     fprintf('Figure %d shows model output with estimated thetaMAP Idnetity S, as well as model output error.\n\n', figure_num);
%     figure(figure_num);
%     figure_num = figure_num + 1;
%     fig = gcf;
%     fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
%     subplot(2,1,1);
%     plot(p3_y);
%     hold on;
%     plot(y_estimated);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'MAP Actual vs. Model Output with S identity';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     lgd = legend('$y(k)$', '$y_{\hat{\theta}_{ML}}(k)$');
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = 12;
%     lgd.Location = 'northwest';
%     subplot(2,1,2);
%     plot(estimation_error);
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'MAP Model Output Error  with S identity';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     axes.YLabel.Interpreter = 'latex';
%     axes.YLabel.String = '$y(k) - y_{\hat{\theta}_{ML}}(k)$';
%     axes.YLabel.FontSize = 14;
% 
%     estimation_error_var = var(estimation_error(:), 1);
%     fprintf('MAP, S identity, Estimated Model ouput error variance is: %.4f.\n', estimation_error_var);
%     fprintf('MAP, S identity, Estimated Model ouput error is: %.4f.\n', mean(estimation_error(:).^2));
    %-----------------------------------------------------------------------
    %Cross Validation
    %-----------------------------------------------------------------------
    Phi = zeros(N, num_params);
    Phi(2,:) = [p3_u_cv(1) 0 0 0 0 0 0 0];
    Phi(3,:) = [p3_u_cv(2) p3_u_cv(1) 0 0 0 0 0 0];
    Phi(4,:) = [p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0 0 0 0];
    Phi(5,:) = [p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0 0 0];
    Phi(6,:) = [p3_u_cv(5) p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0 0];
    Phi(7,:) = [p3_u_cv(6) p3_u_cv(5) p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0 0];
    Phi(8,:) = [p3_u_cv(7) p3_u_cv(6) p3_u_cv(5) p3_u_cv(4) p3_u_cv(3) p3_u_cv(2) p3_u_cv(1) 0];
    for k = num_params+1 : N
        Phi(k,:) = [p3_u_cv(k-1) p3_u_cv(k-2) p3_u_cv(k-3) p3_u_cv(k-4) p3_u_cv(k-5) p3_u_cv(k-6) p3_u_cv(k-7) p3_u_cv(k-8)];
    end
    
    cv_errors = zeros(5, 1);
    
%     figure(figure_num);
%     figure_num = figure_num + 1;
%     fig = gcf;
%     fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
%     subplot(2,1,1);
%     plot(p3_y_cv);
%     hold on;
    for i = 1:5
        y_estimated = Phi * p3_b_MAP(:, i);
        estimation_error = p3_y_cv - y_estimated;
        cv_errors(i) = mean(estimation_error(:).^2);
%         subplot(2,1,1);
%         plot(y_estimated);
%         hold on;

%     subplot(2,1,2);
%     plot(estimation_error);

        
    end
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'MAP Output with different S';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     lgd = legend('$Y(k)$', '$S_{1}(k)$', '$S_{2}(k)$', '$S_{3}(k)$', '$S_{4}(k)$', '$S_{5}(k)$');
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = 12;
%     lgd.Location = 'northwest';
%     for i = 1:5
%         y_estimated = Phi * p3_b_MAP(:, i);
%         estimation_error = p3_y_cv - y_estimated;
%         cv_errors(i) = mean(estimation_error(:).^2);
%         subplot(2,1,2);
%         plot(estimation_error);
%         hold on;
%     end
%     axis tight;
%     axes = gca;
%     axes.Title.Interpreter = 'latex';
%     axes.Title.String = 'MAP Model Error with different S';
%     axes.Title.FontSize = 18;
%     axes.XAxis.TickLabelInterpreter = 'latex';
%     axes.XAxis.FontSize = 10;
%     axes.YAxis.TickLabelInterpreter = 'latex';
%     axes.YAxis.FontSize = 10;
%     axes.XLabel.Interpreter = 'latex';
%     axes.XLabel.String = '$k$ $[samples]$';
%     axes.XLabel.FontSize = 14;
%     lgd = legend('$S_{1}(k)$', '$S_{2}(k)$', '$S_{3}(k)$', '$S_{4}(k)$', '$S_{5}(k)$');
%     lgd.Interpreter = 'latex';
%     lgd.FontSize = 12;
%     lgd.Location = 'northwest';
    
    fprintf('Figure %d shows cross validation to choose the appropiate covariance matrix for the MAP estimate.\n\n', figure_num);
    figure(figure_num);
    figure_num = figure_num + 1;
    fig = gcf;
    fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
    plot(cv_errors);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Cross Validation Error Results';
    axes.Title.FontSize = 18;
    axes.XAxis.TickLabelInterpreter = 'latex';
    axes.XAxis.FontSize = 10;
    axes.YAxis.TickLabelInterpreter = 'latex';
    axes.YAxis.FontSize = 10;
    axes.XLabel.Interpreter = 'latex';
    axes.XLabel.String = '$S_i$';
    axes.XLabel.FontSize = 14;
    axes.YLabel.Interpreter = 'latex';
    axes.YLabel.String = '$(y_{cv}(k) - y_{\theta}(k))^2$';
    axes.YLabel.FontSize = 14;
    
    fprintf('Figure %d shows the effect of different priors on the parameter values.\n\n', figure_num);
    figure(figure_num);
    figure_num = figure_num + 1;
    fig = gcf;
    fig.Position = [mod(figure_num,2)*screenwidth/2, 0, screenwidth/2, screenheight];
    for i = 1:5
        plot(p3_b_MAP(:, i));
        hold on;
    end
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Effect of Priors on Parameter Values';
    axes.Title.FontSize = 18;
    axes.XAxis.TickLabelInterpreter = 'latex';
    axes.XAxis.FontSize = 10;
    axes.YAxis.TickLabelInterpreter = 'latex';
    axes.YAxis.FontSize = 10;
    axes.XLabel.Interpreter = 'latex';
    axes.XLabel.String = '$b_i$';
    axes.XLabel.FontSize = 14;
    axes.YLabel.Interpreter = 'latex';
    axes.YLabel.String = '$\theta$';
    axes.YLabel.FontSize = 14;
    lgd = legend('$S_{1}(k)$', '$S_{2}(k)$', '$S_{3}(k)$', '$S_{4}(k)$', '$S_{5}(k)$');
    lgd.Interpreter = 'latex';
    lgd.FontSize = 12;
    lgd.Location = 'northwest';
    
    p3_cv_error    	= cv_errors;   % vector of dimension 5x1
    [~,Index] = min(p3_cv_error);
    p3_prior_best = Index;            % scalar integer in the set {1,2,3,4,5}


end
