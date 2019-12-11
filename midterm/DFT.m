function UN = DFT(u)

    N = size(u, 2);
    j = sqrt(-1);
    w_n = 0:N-1;
    w_n = w_n * (2*pi/N);
    w_n = repmat(w_n, N, 1);
    k = 0:N-1;
    k = repmat(k', 1, N);
    exponential = exp((-j .* w_n) .* k);
    UN = u * exponential;
end