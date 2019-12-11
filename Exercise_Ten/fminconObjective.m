function f = fminconObjective(x, num_parameters, mu_v)
    v = x(num_parameters+1:end);
    f = (v - mu_v)' * (v - mu_v);
end