function [sinograms] = apply_butterworth_window_to_sinogram(sinograms, n_butter, low_cutoff, high_cutoff)
% APPLY_BUTTERWORTH_WINDOW_TO_SINOGRAM Multiply sinogram with frequency response
% of a butterworth bandpass filter to obtain smooth decays at the sinogram
% boundaries. The magnitude at the cutoffs is 1/sqrt(2).
[z,p,k] = butter(n_butter ,[low_cutoff high_cutoff]./size(sinograms, 1), 'bandpass');
[b,a] = zp2tf(z, p, k);
h=freqz(b, a, size(sinograms,1));
sin_mean = mean(sinograms,1);
sinograms = sin_mean + (sinograms-sin_mean) .* abs(h);
end
