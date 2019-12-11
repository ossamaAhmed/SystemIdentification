function [thetam1, theta_map] = HS2019_SysID_Exercise_10_18936872()
    
    % Initialization
    clc;
    clear all;
    close all;

    % Screen size used to place plots
    screensize = get(groot, 'ScreenSize');
    screenwidth = screensize(3);
    screenheight = screensize(4);
    %load the data
    data = load('SysID_Exercise_10.mat');
    u1 = data.u1;
    u2 = data.u2;
    y1 = data.y1;
    y2 = data.y2;
    % Part 1
    N = length(y1);
    beta1 = 0.5;
    beta2 = 0.2;
    mu_v = 1.5;
    var_v = 1.2;    
    %get the MLE of the data
    % Nonlinear program solver fmincon
    % x = [theta; v]
    num_parameters = 1;
    fun = @(x)fminconObjective(x, num_parameters, mu_v);
    x0 = zeros(num_parameters + N, 1);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    nonlcon = @(x)fminconConstraint(x, num_parameters, y1, u1, beta1, beta2);
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000);
    x1 = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    theta1 = x1(1)
    v1 = x1(num_parameters+1:end);
    y1_est = y1 - v1;
    %plot the error
    figure(1);
    % y Actual vs. Estimate
    subplot(2,1,1);
    plot(y1);
    hold on;
    plot(y1_est);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'y Actual vs. Estimate';
    axes.Title.FontSize = 18;
    % Error
    subplot(2,1,2);
    plot(v1);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Error';
    axes.Title.FontSize = 18;
    %----------------------------------------------------------------------
    %----------------------PART TWO----------------------------------------
    %----------------------------------------------------------------------
    
%     %calculate w
%     
    w = zeros(N, 1);
    w(2) = u2(1);
    w(3) = u2(2) + (2*u2(1));
    for k = 4 : N
        w(k) = u2(k-1) + (2*u2(k-2)) + u2(k-3);
    end
    %closed form solution of theta
    theta2 = (((1/sqrt(var_v)) * ((y2 - mu_v)' * w)) + 6) / (((1/sqrt(var_v)) * (w' * w)) + 4);
    
    %different solution Ax = b
    N = length(y2);
    mu_v = 1.5;
    var_v = 1.2;
    mu_theta_prior = 1.5;
    var_theta_prior = 0.5 ^ 2;
    n_parameters = 1;
    Aeq = zeros(N, n_parameters);
    beq = y2 - mu_v; 
    Aeq(2,:) = u2(1);
    beq(2) = beq(2) - Aeq(2,:) * mu_theta_prior;
    Aeq(3,:) = u2(2) + 2 * u2(1);
    beq(3) = beq(3) - Aeq(3,:) * mu_theta_prior;
    for k = 4 : N
        Aeq(k,:) = u2(k-1) + 2 * u2(k-2) + u2(k-3);
        beq(k) = beq(k) - Aeq(k,:) * mu_theta_prior;
    end
    Aeq = [Aeq eye(N)];
    
    % x = [theta - mTheta; v - mu_v]
    H = [1 / (2 * var_theta_prior) zeros(n_parameters, N); zeros(N, n_parameters) 1 / (2 * var_v) * eye(N)];
    fT = [];
    A = [];
    b = [];
    x2 = quadprog(H, fT, A, b, Aeq, beq);
    theta2 = x2(1) + mu_theta_prior
    v2 = x2(n_parameters+1:end) + mu_v;
    y2_est = y2 - v2;
    
    %plotting time
    %plot the error
    figure(2);
    % y Actual vs. Estimate
    subplot(2,1,1);
    plot(y2);
    hold on;
    plot(y2_est);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'y2 Actual vs. Estimate';
    axes.Title.FontSize = 18;
    % Error
    subplot(2,1,2);
    plot(v2);
    axis tight;
    axes = gca;
    axes.Title.Interpreter = 'latex';
    axes.Title.String = 'Error';
    axes.Title.FontSize = 18;
    thetam1 = theta1;
    theta_map = theta2;
    function f = fminconObjective(x, num_parameters, mu_v)
        v = x(num_parameters+1:end);
        f = (v - mu_v)' * (v - mu_v);
    end

    function [c, ceq] = fminconConstraint(x, num_parameters, y, u, beta1, beta2)
        c = [];
        theta = x(1:num_parameters);
        v = x(num_parameters+1:end);
        M = length(v);
        v_est = zeros(M, 1);
        v_est(1) = y(1);
        v_est(2) = y(2) - theta * u(1) - beta1 * v_est(1);
        for j = 3 : M
            v_est(j) = y(j) - theta * u(j-1) - beta1 * v_est(j-1) - beta2 * v_est(j-2);
        end
        ceq = v - v_est;
    end
end