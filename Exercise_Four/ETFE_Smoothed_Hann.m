function [omega, G_result] = ETFE_Smoothed_Hann(u, y, gamma)
    lags = -gamma:gamma;
    spectral_input = autocorrelation_finite(u, lags);
    cross_spectral = crosscorrelation_finite(y, u, lags);
    N = size(lags, 2);
    [~, wHann] = WtHann(gamma, N);
    windowed_spectral_input = wHann' .* spectral_input;
    windowed_cross_spectral = wHann' .* cross_spectral;
    windowed_spectral_input_fft = fft(windowed_spectral_input);
    windowed_cross_spectral_fft = fft(windowed_cross_spectral);
    G_result = windowed_cross_spectral_fft ./ windowed_spectral_input_fft;
    N_size = size(G_result, 2);
    omega = (2*pi/N_size)*[0:N_size-1];
end


% function Gest = WHtEst(u, y, gamma)
%     t = -gamma : 1 : gamma;
%     Ru = autocorrelation_F(u, t);
%     Ryu = crosscorrelation_F(y, u, t);
%     W = WHtdom(gamma, t);
%     Ruw = W .* Ru;
%     Ryuw = W .* Ryu;
%     Su = DFT(Ruw);
%     Syu = DFT(Ryuw);
%     Gest = Syu ./ Su;
% end