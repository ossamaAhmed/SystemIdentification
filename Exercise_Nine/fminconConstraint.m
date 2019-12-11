function [c, ceq] = fminconConstraint(x, y, u)
    c = [];
    theta = x(1:4);
    %calculate w
    N = length(y);
    w = zeros(N, 1);
    w(2) = u(1)*theta(3);
    phi = zeros(N, 3);
    phi(2, :) = [-w(1) 0 u(1)];
    for k = 3 : N
        w(k) = -w(k - 1)*theta(1) - w(k - 2)*theta(2) + u(k-1)*theta(3);
        phi(k, :) = [-w(k-1) -w(k-2) u(k-1)];
    end
    
    e = x(5:end);
    N_error = length(e);
    phiTe = zeros(N_error, 1);
    for k = 2 : N_error
        phiTe(k) = e(k-1);
    end
    phiT = [phi phiTe];
    ceq = y - phiT * theta - e;
end