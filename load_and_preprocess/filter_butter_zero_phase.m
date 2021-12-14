function [signal_matrix_filt] = filter_butter_zero_phase(signal_matrix,sampling_frequency,filter_frequencies,use_fft)

%% Zero-phase butterworth filtering in f-domain with circular convolution or with filtfilt function
% NOTE:
%  - Utilizing filtering with fft is faster than 'filtfilt'.
%    However, fft-filtering introduces transients at start and end of the 
%    signal, if bad windowing is applied to signals. This can be avoided using
%    'filtfilt' instead.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter in f-domain if not specified differently
if(nargin < 4)
    use_fft = true;
end

% Demean signal (standalone butter does not demean signal appropriately)
if (filter_frequencies(1) > 0)
    signal_matrix_filt = signal_matrix - mean(signal_matrix,1); 
else
    signal_matrix_filt = signal_matrix;
end

% Create butterworth filter
zeros_filter = [];
poles_filter = [];
gain_filter = 1;
if(filter_frequencies(2) > 0) % 8th order butterworth LP filter
    [z_lower,p_lower,k_lower] = butter(8,filter_frequencies(2).*2./sampling_frequency, 'low');       % Doubled frequency, as "butter" requires cutt-off freq. normed by Nyquist rate
    zeros_filter = [zeros_filter; z_lower];
    poles_filter = [poles_filter; p_lower];
    gain_filter = gain_filter*k_lower;
end
if(filter_frequencies(1) > 0) % 2nd order butterworth HP filter
    [z_higher,p_higher,k_higher] = butter(2,filter_frequencies(1).*2./sampling_frequency, 'high');      % Doubled frequency, as "butter" requires cutt-off freq. normed by Nyquist rate
    zeros_filter = [zeros_filter; z_higher];
    poles_filter = [poles_filter; p_higher];
    gain_filter = gain_filter*k_higher;
end
sos_filter = zp2sos(zeros_filter,poles_filter,gain_filter);
filt = dfilt.df2sos(sos_filter);

if use_fft
    % Create filter transfer function
    IR_filter_hat = freqz(filt, size(signal_matrix,1),'whole',sampling_frequency); % butterworth filter
    IR_filter_hat_zero_phase = abs(IR_filter_hat).^2; % zero-phase butterworth filter
    
    % Filter signal along 1st axis in f-domain
    signal_matrix_filt = ifft(IR_filter_hat_zero_phase.*fft(signal_matrix_filt),'symmetric');
else
    % Filter signal along 1st axis with zero-phase filter function
    signal_matrix_filt = filtfilt(filt.sosMatrix,filt.ScaleValues,signal_matrix_filt);
end
end

