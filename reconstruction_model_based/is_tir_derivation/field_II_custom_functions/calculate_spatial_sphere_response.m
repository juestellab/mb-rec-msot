function [signals_nShape_sir,travel_times] = calculate_spatial_sphere_response(model,filt,filter_length,frequency_sampling_field,detect_single_aperture,sphere_diam,phantom_loc,transducer_number)
%% Compute OA reference signal convolved with SIR for all detectors
%% TODO: Update 3D refraction scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% signal matrix from reference N-shape convolved with computed SIR
if isempty(transducer_number)
    signals_nShape_sir = zeros(model.Data.numTimeSamplesPerTransducer,model.Probe.detector.numOfTransducers);
    transducer_number = 1:model.Probe.detector.numOfTransducers;
else
    signals_nShape_sir = zeros(model.Data.numTimeSamplesPerTransducer,1);
end
    
% time variables
samples_missing_prior_signal=model.Probe.DAC.delayBeforeRecording+model.DataPreprocessing.numCroppedSamplesAtSinogramStart;
time_vec_dac = double((0:(model.Data.numTimeSamplesPerTransducer-1))+samples_missing_prior_signal)'./model.Probe.DAC.frequency;                 

% Virtual source location of phantom for all detectors
if model.useSingleSpeedOfSound
    speed_of_sound_for_SIR = model.Medium.speedOfSound;
     
    x_phantom_virtuals = repmat(phantom_loc(1),[length(transducer_number),1]);
    y_phantom_virtuals = repmat(phantom_loc(2),[length(transducer_number),1]);
      
    travel_times = sqrt(x_phantom_virtuals.^2+y_phantom_virtuals.^2)./model.Medium.speedOfSound;
else
    speed_of_sound_for_SIR = model.Probe.coupling.speedOfSound;
       
    x_transducers = model.Probe.detector.xPositionsOfTransducers(transducer_number);
    y_transducers = model.Probe.detector.yPositionsOfTransducers(transducer_number);
    [x_membrane_intersections,travel_times,~] = refraction(model.Probe.membrane,model.Probe.detector.radius,phantom_loc(1),phantom_loc(2),...
        x_transducers,y_transducers, model.Medium.speedOfSound,speed_of_sound_for_SIR);
    y_membrane_intersections = model.Probe.membrane(x_membrane_intersections);
    x_phantom_virtuals = x_transducers + (x_membrane_intersections - x_transducers)./sqrt((x_membrane_intersections - x_transducers).^2 + (y_membrane_intersections - y_transducers).^2).*travel_times*speed_of_sound_for_SIR;
    y_phantom_virtuals = y_transducers + (y_membrane_intersections - y_transducers)./sqrt((x_membrane_intersections - x_transducers).^2 + (y_membrane_intersections - y_transducers).^2).*travel_times*speed_of_sound_for_SIR;
    nan_indices = isnan(x_membrane_intersections(:,:,1));
               
    % If xmin is not defined, set the virtual location to the actual location of the pixel
    x_phantom_virtuals(isnan(x_membrane_intersections)) = repmat(model.Discretization.region.xPositionsFov(nan_indices),[length(x_transducers),1]);
    y_phantom_virtuals(isnan(x_membrane_intersections)) = repmat(model.Discretization.region.yPositionsFov(nan_indices),[length(y_transducers),1]);
end
    
% N shaped OA signal (reference signal for spherical source of corresponding size)
length_n_shape = max(filter_length,ceil(sphere_diam/speed_of_sound_for_SIR*frequency_sampling_field)+1);
t_vec_nShape = (-floor(length_n_shape/2):floor(length_n_shape/2))/frequency_sampling_field;
nShape = -(t_vec_nShape*speed_of_sound_for_SIR).*((abs(2*t_vec_nShape*speed_of_sound_for_SIR)-sphere_diam) < 1e-9);
nShape = nShape'/sphere_diam*2;                    % Rescaled for normalized model
    
% Iterate over detectors
for index_transducer = transducer_number
    	
    % Homogeneous transformation to have source position wrt to detector
    transducer_angle = ((180-model.Probe.detector.angularCoverage)/2 + (index_transducer-1)*model.Probe.detector.angularCoverage/(model.Probe.detector.numOfTransducers-1))*pi/180;
    theta_rot=transducer_angle+pi/2;
    H_mat = [cos(-theta_rot) -sin(-theta_rot) 0;
             sin(-theta_rot) cos(-theta_rot)  model.Probe.detector.radius;
             0 0 1];
    if length(transducer_number) > 1
        pos_virtual_rot = H_mat*[x_phantom_virtuals(index_transducer);y_phantom_virtuals(index_transducer);1];
    else
        pos_virtual_rot = H_mat*[x_phantom_virtuals;y_phantom_virtuals;1];
    end
      
    % SIR & travel time with FIELD II
    [h_SIR,t_start_SIR] = calc_h(detect_single_aperture,[pos_virtual_rot(1) phantom_loc(3) pos_virtual_rot(2)]);
    h_SIR(h_SIR<1) = 0;                                   % set all small values to 0
        
    % Convolve with n-shape (& apply anti-aliasing filter)
    %% TODO: Does filtering work as expected
    fft_length = length(h_SIR) + length_n_shape - 1;
    n_shape_hat = fft(nShape,fft_length);
    h_f_hat = abs(freqz(filt,fft_length,'whole',frequency_sampling_field)).^2;
    h_SPR_filt = ifft(fft(h_SIR,fft_length,1).*n_shape_hat.*h_f_hat,fft_length,1,'symmetric');
    
    % resample SPR with numerical derivative of antiderivative
    time_vec_field =  t_start_SIR + (-floor(length_n_shape/2):(length(h_SIR)-1+floor(length_n_shape/2)))'/frequency_sampling_field;
    H_SPR_filt = cumtrapz(time_vec_field,h_SPR_filt);
    for index_time_vec_dac = 1:length(time_vec_dac)
        [~,index_time_vec_field] = min(abs(time_vec_field-time_vec_dac(index_time_vec_dac)));
        if (index_time_vec_field <= length(time_vec_field))
            if (time_vec_dac(index_time_vec_dac) <= time_vec_field(index_time_vec_field))
                if index_time_vec_field == 1
                    if 1 < length(transducer_number)
                        signals_nShape_sir(index_time_vec_dac,index_transducer) = H_SPR_filt(index_time_vec_field)*frequency_sampling_field;
                    else
                        signals_nShape_sir(index_time_vec_dac) = H_SPR_filt(index_time_vec_field)*frequency_sampling_field;
                    end
                else
                    if 1 < length(transducer_number)
                        signals_nShape_sir(index_time_vec_dac,index_transducer) = (H_SPR_filt(index_time_vec_field)-H_SPR_filt(index_time_vec_field-1))*frequency_sampling_field;
                    else
                        signals_nShape_sir(index_time_vec_dac) = (H_SPR_filt(index_time_vec_field)-H_SPR_filt(index_time_vec_field-1))*frequency_sampling_field;
                    end
                end
            else
                if (index_time_vec_field < length(time_vec_field))
                    if 1 < length(transducer_number)
                        signals_nShape_sir(index_time_vec_dac,index_transducer) = (H_SPR_filt(index_time_vec_field+1)-H_SPR_filt(index_time_vec_field))*frequency_sampling_field;
                    else
                        signals_nShape_sir(index_time_vec_dac) = (H_SPR_filt(index_time_vec_field+1)-H_SPR_filt(index_time_vec_field))*frequency_sampling_field;
                    end
                end
            end
        end
    end
end
    
end