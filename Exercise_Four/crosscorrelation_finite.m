function result = crosscorrelation_finite(y, x, shifts)
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
    R_y_x = y * x_shifted;
    result = R_y_x;
end