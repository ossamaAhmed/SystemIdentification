function result = autocorrelation_periodic(x, shifts)
    N = size(x, 2); %period
    num_shifts = size(shifts, 2);
    x_shifted = zeros(N, num_shifts);
    for i=1:num_shifts
        %shift the period signal
        x_shifted(:, i) = circshift(x, shifts(1, i))';
    end
    R_x = 1/N * x * x_shifted;
    result = R_x;
    
end