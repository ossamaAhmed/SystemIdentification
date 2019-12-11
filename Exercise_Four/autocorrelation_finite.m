function result = autocorrelation_finite(x, shifts)
     N = size(x, 2); %period
    num_shifts = size(shifts, 2);
    x_shifted = zeros(N, num_shifts);
    for i=1:num_shifts
        %shift the period signal
        s = shifts(1,i);
        if s >= 0
            x_shifted(s+1:N,i) = x(1,1:N-s)';
        else
            x_shifted(1:N+s,i) = x(1,-s+1:N);
        end
    end
    R_x = x * x_shifted;
    result = R_x;
end
