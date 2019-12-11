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


% function [G_result] = ETFE_Smoothed_Hann(u, y, gamma)
%     U_freq = fft(u);
%     Y_freq = fft(y);
%     G_estimated = Y_freq ./ U_freq;
%     N = size(G_estimated, 2);
%     [omega, Wg] = WfHann(gamma, N);
%     zidx = find(Wg==0);
%     omega = [omega(zidx:N); omega(1:zidx-1)];
%     Wg = [Wg(zidx:N), Wg(1:zidx-1)];
%     a = U_freq.*conj(U_freq);
%     G_result = zeros(1, N);
%     for wn =1:N
%         Wnorm = 0;
%         for xi = 1:N
%             widx = mod(xi-wn, N)+1;
%             G_result(wn) = G_result(wn) + (Wg(widx) .* G_estimated(xi) .* a(xi));
%             Wnorm = Wnorm + (Wg(widx) .* a(xi));
%         end
%         G_result(wn) = G_result(wn) ./ Wnorm;
%     end
% end

% 
% function Gn = GestWHf(u, y, gamma)
%     U = DFT(u);
%     Uabs2 = abs(U) .^ 2;
%     Y = DFT(y);
%     Gn1 = Y ./ U;
%     m = size(Gn1, 2);
%     n = size(WHf, 2);
%     Gn = zeros(1,m);
%     for i = 1:m
%         Wsum = 0;
%         for j = 1:n
%             Gn(1,i) = Gn(1,i) + WHf(1,j) * Uabs2(1,i-j) * Gn1(1,i-j);
%             Wsum = Wsum + WHf(1,j) * Uabs2(1,i-j);
%         end
% %     end
% end


% function [G_result] = ETFE_Smoothed_Hann(u, y, gamma)
%     U_freq = fft(u);
%     Y_freq = fft(y);
%     G_estimated = Y_freq ./ U_freq;
%     N = size(G_estimated, 2);
%     [~, Wg] = WfHann(gamma, N);
%     a = 1/N .* (abs(U_freq) .^2);
%     weighted_G_estimate = G_estimated .* a;
%     G_result = conv(weighted_G_estimate, Wg, 'same');
%     norm = conv(a, Wg, 'same');
%     G_result = G_result ./ norm;
% end

%     G1_filt_num = zeros(period_len,num_periods);
%     G1_filt_den = zeros(period_len,num_periods);
% 
%     % Filter each period of the ETFE estimate, using frequency-domain
%     % windowing as on slide 4.17
%     for freq = 1:period_len
%         modified_freq_inds = mod((1:period_len) - freq,period_len)+1;
% 
%         % Calculate discretized integrals
%         for period_num = 1:num_periods
%             G1_filt_num(freq,period_num) = sum(freq_window(modified_freq_inds).*abs(U1).^2.*G1_mat(:,period_num));
%             G1_filt_den(freq,period_num) = sum(freq_window(modified_freq_inds).*abs(U1).^2);
%         end
%     end
%     
%    
%     
%     function Gn = GestWHf(u, y, gamma)
%     U = DFT(u);
%     Uabs2 = abs(U) .^ 2;
%     Y = DFT(y);
%     Gn1 = Y ./ U;
%     m = size(Gn1, 2);
%     n = size(WHf, 2);
%     Gn = zeros(1,m);
%     for i = 1:m
%         Wsum = 0;
%         for j = 1:n
%             Gn(1,i) = Gn(1,i) + WHf(1,j) * Uabs2(1,i-j) * Gn1(1,i-j);
%             Wsum = Wsum + WHf(1,j) * Uabs2(1,i-j);
%         end
% %     end
% end

    % Compute filtered transfer function by averaging each period's filtered ETFE
%     G1_filt = G1_filt_num./G1_filt_den;