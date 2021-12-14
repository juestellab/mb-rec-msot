function [data_crop] = filter_and_crop_MSOT_data(data_raw, daq_frequency, first_relevant_signal_per_detector, last_relevant_signal_per_detector, filt_min, filt_max, skipped_prefix)
%% Filter and Crop MSOT measurement data
% Input: data_raw - data loaded from MSOT file with layout [samples, detectors, wavelengths, frames, slices, runs]
%        daq_frequency - sampling frequency of DAQ in probe
%        first_relevant_signal_per_detector
%        last_relevant_signal_per_detector
%        filt_min - low cutoff frequency of filter
%        filt_max - high cutoff frequency of filter
%        skipped_prefix - number of samples that are already removed from the beginning of data_raw
%
% Output: data_filt_crop - filtered and cropped MSOT measurement data
%% set defaults
if (nargin < 7) || isempty(skipped_prefix)
    skipped_prefix = 0;
end
if (nargin < 6) || isempty(filt_max)
    filt_max = 8e6;
end
if (nargin < 5) || isempty(filt_min)
    filt_min = 0.05e6;
end

% Apply butterworth BP filter
[data_filt] = filter_butter_zero_phase(data_raw, daq_frequency, [filt_min, filt_max]);

% Crop SIR to adequate length
[data_crop] = crop_MSOT_data(data_filt, first_relevant_signal_per_detector, last_relevant_signal_per_detector, skipped_prefix);
