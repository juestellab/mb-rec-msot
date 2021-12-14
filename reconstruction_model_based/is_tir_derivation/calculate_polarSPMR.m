function [SPMR_polar,SPMR_timestep_index_start_and_end]  = calculate_polarSPMR(probe_id,spmr_version,speed_of_sound_coupling, transducer_elevation_radius, transducer_height, transducer_pitch, inter_transducer_spacing, dac_sampling_frequency, distance_fov_min, distance_fov_max, distance_fov_stepsize, angle_azimuth_fov_max, angle_azimuth_fov_stepsize, elevation_fov_max, elevation_fov_stepsize, pixel_model_function, pixel_model_hash)
                              
%% Calculation of polar Spatial Pixel Model Response lookup table for a single detector element
% NOTES: 
%  - FIELD II SIR computation is based on following paper:
%    J.A. Jensen, N.B. Svendsen. "Calculation of Pressure Fields from Arbitrarily Shaped, Apodized, and Excited Ultrasound Transducers." 
%    IEEE Transactions on Ultrasonics, Ferroelectrics and Frequency Control, vol. 39, no. 2, pp. 262-267, 1991
%  - FIELD II does not allow usage of parallelization toolbox ("parfor")
%  - Outer bounds of polar grid can be estimated by drawing a circle through edges of field of view which has 
%    to be fully covered by polar field of view.
%    calculate_SIR_params() is a method which can be used to estimate the bounds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check, if 2D or 3D is computed
if (elevation_fov_max > 0)
    fprintf('calculate_polarSPMR: 3D polar SPMR for single transducer of detector array of probe %s will be computed.\n', probe_id); 
else
    fprintf('calculate_polarSPMR: 2D polar SPMR for single transducer of detector array of probe %s will be computed.\n', probe_id); 
end

% Compute polar grid parameters
distance_fov_end_index = ceil((distance_fov_max-distance_fov_min)/distance_fov_stepsize);
distance_mat_SPMR = distance_fov_min+(0:distance_fov_end_index)*distance_fov_stepsize;   % radial distance (for IBMI coord. frame: r=x.^2+z.^2)
angle_azimuth_fov_end_index = ceil(angle_azimuth_fov_max/angle_azimuth_fov_stepsize);
angle_azimuth_mat_SPMR = (0:angle_azimuth_fov_end_index)*angle_azimuth_fov_stepsize;     % azimuth angle (for IBMI coord. frame: theta = atan(x/z))
elevation_fov_end_index = ceil(elevation_fov_max/elevation_fov_stepsize);
elevation_mat_SPMR_full = (0:elevation_fov_end_index)*elevation_fov_stepsize;            % out-of-focus plane elevation (for IBMI coord frame: y)

% Find starting index of polarSIR on elevation axis & limit maximum number of elevation slices per polarSIR to 10
elevation_mat_length_max = 10;
elevation_mat_start_index = 1;
elevation_mat_end_index = min(length(elevation_mat_SPMR_full),elevation_mat_start_index+elevation_mat_length_max-1);
                    
% Init FIELD II, set params & define transducer model
field_sampling_frequency = 5e9;
elevation_focused_transducer = init_field_II_and_get_single_toroidal_focused_transducer(transducer_elevation_radius,transducer_height,transducer_pitch,inter_transducer_spacing,distance_fov_min,speed_of_sound_coupling,field_sampling_frequency);

% Files & paths
path_to_current_folder = fileparts(which('calculate_polarSPMR'));     % global path
filename = strcat('polarSPMR_',probe_id,'_c_',num2str(speed_of_sound_coupling),'_pixelModel_',num2str(pixel_model_hash));

% Time variable index bounds
% Note: Length of SIR is fixed initially to compute the IR of the filter only once (Filter-IR in every iteration is very inefficient)
half_number_time_samples_SPMR_dac_sampling = 45;        % This value is usually sufficient for a 2D SIR and the field of view, we are using    
upsampling_factor = ceil(field_sampling_frequency/dac_sampling_frequency);
half_number_time_samples_SIR_field_sampling = upsampling_factor*half_number_time_samples_SPMR_dac_sampling; % Less samples can be used for padding the FIELD signal, as the output is not centered around the travel time
           
% filtering type
filter_type = 'Butterworth';
filter_cutoff = dac_sampling_frequency/2;
filter_order = 8;
[zeros_AntiAlias,poles_AntiAlias,gain_AntiAlias] = butter(8,filter_cutoff.*2./field_sampling_frequency, 'low');      % doubled frequency, as "butter" requires cut-off freq. normed by Nyquist rate
sos_filter = zp2sos(zeros_AntiAlias,poles_AntiAlias,gain_AntiAlias);
filter_AntiAlias = dfilt.df2sos(sos_filter);

% N-shape pixel model solution of wave equation
% (ensure that pixel model decays to 0 at some point)
half_number_time_samples_filter_field_sampling = 200;
if ~(abs(pixel_model_function(half_number_time_samples_filter_field_sampling/field_sampling_frequency*speed_of_sound_coupling)) < 1e-3)
    time_sample_bounds = (1:0.25:5)*1e-7;
    for time_sample_bound = time_sample_bounds
        if (abs(pixel_model_function(time_sample_bound*speed_of_sound_coupling)) < 1e-3)
            half_number_time_samples_filter_field_sampling = ceil(time_sample_bound*field_sampling_frequency);
            break;
        end
    end
end
time_samples_filt = (-half_number_time_samples_filter_field_sampling:half_number_time_samples_filter_field_sampling)'/field_sampling_frequency;
pixel_model = (pixel_model_function(time_samples_filt*speed_of_sound_coupling).*(-time_samples_filt));
pixel_model_string_SPMR = func2str(pixel_model_function);

% compute fft of signals
fft_length = 2*(half_number_time_samples_SIR_field_sampling+half_number_time_samples_filter_field_sampling)-1;
IR_filter_hat = abs(freqz(filter_AntiAlias,int32(fft_length),'whole',field_sampling_frequency)).^2; 
pixel_model_hat = fft(pixel_model,int32(fft_length),1);

% resampling time vectors
time_samples_dac = (-half_number_time_samples_SPMR_dac_sampling:half_number_time_samples_SPMR_dac_sampling)'/dac_sampling_frequency;
time_samples_field =  ((1:fft_length)-half_number_time_samples_filter_field_sampling-1)/field_sampling_frequency;
 
% build the model
SPMR_polar = zeros(min(elevation_mat_end_index-elevation_mat_start_index+1,elevation_mat_length_max),length(angle_azimuth_mat_SPMR),length(distance_mat_SPMR),2*half_number_time_samples_SPMR_dac_sampling+1);
start_end_time_samples_SPMR = cat(4,ones(size(SPMR_polar,1:3)),(2*half_number_time_samples_SPMR_dac_sampling+1)*ones(size(SPMR_polar,1:3)));
computation_time_single_azimuth_angle = 0;
computation_time_single_elevation_slice = 0;
tic;
for index_elevation_mat = elevation_mat_start_index:length(elevation_mat_SPMR_full)
    fprintf('current elevation step (total: %d): %d\n',length(elevation_mat_SPMR_full),index_elevation_mat);
    fprintf('current azimuth angle step (total: %d) (distance steps per angle step: %d): ', length(angle_azimuth_mat_SPMR), length(distance_mat_SPMR));
    
    for index_angle_azimuth_mat = 1:length(angle_azimuth_mat_SPMR)
        if (index_angle_azimuth_mat == 1) || (~mod(index_angle_azimuth_mat,100))
            fprintf('%d, ',index_angle_azimuth_mat);
        end
    
        for index_distance_mat = 1:length(distance_mat_SPMR)
        
            % Get current point
            % Note: Theta is measured from the z-axis (IBMI definition)!
            point = [distance_mat_SPMR(index_distance_mat)*sin(angle_azimuth_mat_SPMR(index_angle_azimuth_mat)*pi/180),...
                     elevation_mat_SPMR_full(index_elevation_mat), ...
                     distance_mat_SPMR(index_distance_mat)*cos(angle_azimuth_mat_SPMR(index_angle_azimuth_mat)*pi/180)];
                    
            % Get travelling time 
            travel_time = norm(point,2)/speed_of_sound_coupling;
            time_samples_resampling = travel_time + time_samples_dac;

            % Compute SIR with FIELD II
            [SIR_point,start_time_SIR_point] = calc_h(elevation_focused_transducer,point);
            SIR_point(SIR_point<1) = 0;                                   % set all small values to 0
            
            % Pad SIR to correct length for convolution & increase size of
            % SPMR_polar if necessary
            zeros_to_be_padded = double(2*half_number_time_samples_SIR_field_sampling+1)-length(SIR_point);
            if zeros_to_be_padded < half_number_time_samples_filter_field_sampling
                fprintf('\n');
                warning('half_number_time_samples_SPMR_dac_sampling has been chosen too short initially and had to be increased during runtime.');
                fprintf('current azimuth angle step (total: %d) (distance steps per angle step: %d): %d, ',length(angle_azimuth_mat_SPMR), length(distance_mat_SPMR), index_angle_azimuth_mat);
                zeros_to_be_padded = zeros_to_be_padded - half_number_time_samples_filter_field_sampling;
                samples_padded_SPMR_dac_sampling = round(abs(zeros_to_be_padded)/upsampling_factor)+5;  % Include buffer of 5 additional samples
                samples_padded_SIR_field_sampling = samples_padded_SPMR_dac_sampling*upsampling_factor;
                half_number_time_samples_SIR_field_sampling = half_number_time_samples_SIR_field_sampling + samples_padded_SIR_field_sampling;
                half_number_time_samples_SPMR_dac_sampling = half_number_time_samples_SPMR_dac_sampling + samples_padded_SPMR_dac_sampling;
                
                % Recompute new fft-length & time sample vectors
                fft_length = 2*(half_number_time_samples_SIR_field_sampling+half_number_time_samples_filter_field_sampling)-1;
                IR_filter_hat = abs(freqz(filter_AntiAlias,int32(fft_length),'whole',field_sampling_frequency)).^2; 
                pixel_model_hat = fft(pixel_model,int32(fft_length),1);
                time_samples_dac = (-half_number_time_samples_SPMR_dac_sampling:half_number_time_samples_SPMR_dac_sampling)'/dac_sampling_frequency;
                time_samples_field =  ((1:fft_length)-half_number_time_samples_filter_field_sampling-1)/field_sampling_frequency;
                time_samples_resampling = travel_time + time_samples_dac;

                % Pad polar_SPMR & adjust starting values
                SPMR_polar = padarray(SPMR_polar,[0,0,0,samples_padded_SPMR_dac_sampling],0,'both'); 
                start_end_time_samples_SPMR = start_end_time_samples_SPMR + samples_padded_SPMR_dac_sampling;
                
                % Compute new padding length of SIR
                zeros_to_be_padded = double(2*half_number_time_samples_SIR_field_sampling+1)-length(SIR_point);
            end
            SIR_point = padarray(SIR_point,[zeros_to_be_padded,0],0,'post'); 
            
            % Convolve with n-shape & apply anti-aliasing filter
            SPMR_filt = ifft(fft(SIR_point,int32(fft_length),1).*pixel_model_hat.*IR_filter_hat,int32(fft_length),1,'symmetric');
            
            % Resample SPMR with 1D-linear interpolation
            time_samples_SPMR = start_time_SIR_point + time_samples_field';
            indices_to_be_sampled = find((time_samples_resampling-time_samples_SPMR(1)) > 0);
            SPMR_filt_dac = interp1(time_samples_SPMR,SPMR_filt,time_samples_resampling(indices_to_be_sampled(1):end),'linear');
            
            % Set small values to 0
            SPMR_filt_dac_max = max(abs(SPMR_filt_dac));
            SPMR_filt_dac(abs(SPMR_filt_dac) < SPMR_filt_dac_max*1e-4) = 0;
            
            % Assign SPMR
            SPMR_polar(index_elevation_mat-elevation_mat_start_index+1,index_angle_azimuth_mat,index_distance_mat,indices_to_be_sampled(1):end)= SPMR_filt_dac;
            
            % Assign start and end indices
            non_zero_indices = find(SPMR_polar(index_elevation_mat-elevation_mat_start_index+1,index_angle_azimuth_mat,index_distance_mat,:));
            start_end_time_samples_SPMR(index_elevation_mat-elevation_mat_start_index+1,index_angle_azimuth_mat,index_distance_mat,1) = non_zero_indices(1);
            start_end_time_samples_SPMR(index_elevation_mat-elevation_mat_start_index+1,index_angle_azimuth_mat,index_distance_mat,2) = non_zero_indices(end);
        end
                
        if (computation_time_single_azimuth_angle == 0)
            computation_time_single_azimuth_angle = toc;
        end
    end
    fprintf('\n');
    
    if (computation_time_single_elevation_slice == 0)
        computation_time_single_elevation_slice = toc;
        fprintf('Computation time for single elevation slice: %d s\n',computation_time_single_elevation_slice);
    end
    
    % save SIR to file
    if index_elevation_mat == elevation_mat_end_index
        fprintf('SPMR will now be saved to a file.\n');
        
        % Crop SPMR to suitable length
        max_half_number_time_samples_SPMR_dac_sampling = max(half_number_time_samples_SPMR_dac_sampling+1-min(vec(start_end_time_samples_SPMR(:,:,:,1))),max(vec(start_end_time_samples_SPMR(:,:,:,2)))-half_number_time_samples_SPMR_dac_sampling-1)+1;       % Add one sample on the bound
        if (max_half_number_time_samples_SPMR_dac_sampling) < half_number_time_samples_SPMR_dac_sampling
            SPMR_polar = SPMR_polar(:,:,:,(-max_half_number_time_samples_SPMR_dac_sampling:max_half_number_time_samples_SPMR_dac_sampling)+half_number_time_samples_SPMR_dac_sampling+1);
            start_end_time_samples_SPMR = start_end_time_samples_SPMR-(half_number_time_samples_SPMR_dac_sampling-max_half_number_time_samples_SPMR_dac_sampling);
        end
        SPMR_timestep_index_start_and_end = int32(start_end_time_samples_SPMR);
            
        if index_elevation_mat == elevation_mat_start_index
            elevation_mat_SPMR = elevation_mat_SPMR_full;
            save(fullfile(path_to_current_folder, '..','data', 'SPMRs', filename),...
                 'SPMR_polar','SPMR_timestep_index_start_and_end','speed_of_sound_coupling','distance_mat_SPMR','angle_azimuth_mat_SPMR','elevation_mat_SPMR','filter_type','filter_cutoff','filter_order','spmr_version','pixel_model_string_SPMR','field_sampling_frequency','-v7.3');
        else
            elevation_mat_SPMR = elevation_mat_SPMR_full(elevation_mat_start_index:elevation_mat_end_index);
            save(fullfile(path_to_current_folder, '..','data', 'SPMRs', strcat(filename,'_elevIndexMin',num2str(elevation_mat_start_index),'_elevIndexMax_',num2str(elevation_mat_end_index))),...
                 'SPMR_polar','SPMR_timestep_index_start_and_end','speed_of_sound_coupling','distance_mat_SPMR','angle_azimuth_mat_SPMR','elevation_mat_SPMR','filter_type','filter_cutoff','filter_order','spmr_version','pixel_model_string_SPMR','field_sampling_frequency','-v7.3');
            elevation_mat_start_index = index_elevation_mat+1;
            elevation_mat_end_index = min(length(elevation_mat_SPMR_full),elevation_mat_start_index+elevation_mat_length_max-1);
            SPMR_polar=zeros(min(elevation_mat_end_index-elevation_mat_start_index+1,elevation_mat_length_max),length(angle_azimuth_mat_SPMR),length(distance_mat_SPMR),2*half_number_time_samples_SPMR_dac_sampling+1);
        end

    end
end

computation_time_SPMR = toc;
fprintf('SPMR computation finished.\n');
fprintf('  Computation time for single azimuth angle: %d s\n',computation_time_single_azimuth_angle);
fprintf('  Computation time for single elevation slice: %d s\n',computation_time_single_elevation_slice);
fprintf('  Total computation time: %d s\n',computation_time_SPMR);            

end
