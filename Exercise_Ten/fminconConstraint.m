function [c, ceq] = fminconConstraint(x, num_parameters, y, u, beta1, beta2)
    c = [];
    theta = x(1:num_parameters);
    v = x(num_parameters+1:end);
    N = length(v);
    v_est = zeros(N, 1);
    v_est(1) = y(1);
    v_est(2) = y(2) - theta * u(1) - beta1 * v_est(1);
    for k = 3 : N
        v_est(k) = y(k) - theta * u(k-1) - beta1 * v_est(k-1) - beta2 * v_est(k-2);
    end
    ceq = v - v_est;
end