classdef MSOTmodel
    properties
        Probe;
        Discretization;
        Medium;
        Data;
        Funcs;
        includeSir;
        useSingleSpeedOfSound;
        DataPreprocessing;
    end
    properties (Constant)
        VALID_SPMR_VERSION = '2021-03-05';
        SPMR_GRID_NUMERICAL_TOLERANCE = 1e-10;
    end
    
    methods (Access=public)
        function model = MSOTmodel(probe, discretization, speed_of_sound_below_membrane, include_sir, use_single_speed_of_sound,num_cropped_samples_at_sinogram_start, filt_cutoff_min, filt_cutoff_max, model_normalization_factor)            
            if nargin < 9
                model_normalization_factor = [];
            end
            
            model.Probe = probe;
            model.Discretization = discretization;
            model.Medium.speedOfSound = speed_of_sound_below_membrane;
            model.includeSir = include_sir; 
            model.useSingleSpeedOfSound = use_single_speed_of_sound;
            model.DataPreprocessing.numCroppedSamplesAtSinogramStart = num_cropped_samples_at_sinogram_start;
            model.DataPreprocessing.filtCutoffMin = filt_cutoff_min;
            model.DataPreprocessing.filtCutoffMax = filt_cutoff_max;
            model.Data.numTimeSamplesPerTransducer = int32(model.Probe.DAC.numRecordedSamplesPerTransducer - model.DataPreprocessing.numCroppedSamplesAtSinogramStart);
            
            [model.Data.absoluteTravelTimeIndices,acoustic_data] = model.calculate_travel_time_indices();
            
            % Calculate radius for circular interpolation during model evaluation
            model.Data.circularInterpolationRadiusInSamples = ...
                (model.Discretization.region.zPositionsFov(2) - model.Discretization.region.zPositionsFov(1))...
                /(model.Medium.speedOfSound/model.Probe.DAC.frequency);

            if model.includeSir
                [model.Data.IRHat,model.Data.IRFlippedHat] = model.get_impulse_response_for_model_with_spmr();
                [model.Data.SPMR,model.Data.SPMRTimestepIndexStartAndEnd,model.Data.SPMRIndices] = model.get_SPMR_and_SPMR_indices(acoustic_data);  
                model.Funcs = model.generate_normalized_model_functions(@model.forward, @model.transpose, model_normalization_factor);
            else
                [model.Data.IRHat,model.Data.IRFlippedHat] = model.get_impulse_response_for_model_without_spmr();
                model.Funcs = model.generate_normalized_model_functions(@model.forwardWithoutSIR, @model.transposeWithoutSIR, model_normalization_factor);
            end
        end
        
        function visualize(model)
            model.Probe.visualize();
            fig = get(groot,'CurrentFigure');

            % FoV displayed
            colorCode_green = [0 0.5 0];
            hold on;
            line([model.Discretization.region.xPositionsFov(1,1) model.Discretization.region.xPositionsFov(1,model.Discretization.sizeOfPixelGrid(1))],...
                 [model.Discretization.region.zPositionsFov(1,1) model.Discretization.region.zPositionsFov(1,model.Discretization.sizeOfPixelGrid(1))], ...
                 'Color',colorCode_green);
            line([model.Discretization.region.xPositionsFov(1,model.Discretization.sizeOfPixelGrid(1)) model.Discretization.region.xPositionsFov(model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1))],...
                 [model.Discretization.region.zPositionsFov(1,model.Discretization.sizeOfPixelGrid(1)) model.Discretization.region.zPositionsFov(model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1))], ...
                 'Color',colorCode_green);
            line([model.Discretization.region.xPositionsFov(model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1)) model.Discretization.region.xPositionsFov(model.Discretization.sizeOfPixelGrid(2),1)],...
                 [model.Discretization.region.zPositionsFov(model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1)) model.Discretization.region.zPositionsFov(model.Discretization.sizeOfPixelGrid(2),1)], ...
                 'Color',colorCode_green);
            line([model.Discretization.region.xPositionsFov(model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1)) model.Discretization.region.xPositionsFov(model.Discretization.sizeOfPixelGrid(2),1)],...
                 [model.Discretization.region.zPositionsFov(model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1)) model.Discretization.region.zPositionsFov(model.Discretization.sizeOfPixelGrid(2),1)], ...
                 'Color',colorCode_green);
            line([model.Discretization.region.xPositionsFov(model.Discretization.sizeOfPixelGrid(2),1) model.Discretization.region.xPositionsFov(1,1)],...
                 [model.Discretization.region.zPositionsFov(model.Discretization.sizeOfPixelGrid(2),1) model.Discretization.region.zPositionsFov(1,1)], ...
                 'Color',colorCode_green);
            
            hold off;
        end
        
        function plotSensitivityMap(model,number_detector,frequency_band_cutoffs)
            if (nargin < 3)
               frequency_band_cutoffs = [model.DataPreprocessing.filtCutoffMin, model.DataPreprocessing.filtCutoffMax]; 
            end
            if (nargin < 2) || isempty(number_detector)
                number_detector = 1:model.Probe.detector.numOfTransducers; 
                fprintf('MSOTModel - plotSensitivityMap(): Sensitivity map of entire detector array is plotted.\n');
            else
                fprintf('MSOTModel - plotSensitivityMap(): Sensitivity map of detector %i is plotted.\n',number_detector);
            end

            % Check if SIR is available & number of detector exists
            assert(~isempty(model.Data.SPMR),'No SPMR included')
            
            % Prepare variables
            TPMR_sensitivity_map = zeros(model.Discretization.numberOfPixels,length(number_detector));
            fft_length = size(model.Data.SPMR,1)+size(model.Probe.detector.EIR,1)-1; 
                        
            % Prepare filter for frequency band
            if ~isempty(frequency_band_cutoffs)
                fprintf('The map is computed for the following frequency band: [%d,%d]\n',max(frequency_band_cutoffs(1),0),min(frequency_band_cutoffs(2),model.Probe.DAC.frequency/2));
                zeros_filter = [];
                poles_filter = [];
                gain_filter = 1;
                if (frequency_band_cutoffs(2) < model.Probe.DAC.frequency/2) % 8th order butterworth LP filter
                    [z_lower,p_lower,k_lower] = butter(8,frequency_band_cutoffs(2).*2./model.Probe.DAC.frequency, 'low');       % Doubled frequency, as "butter" requires cutt-off freq. normed by Nyquist rate
                    zeros_filter = [zeros_filter; z_lower];
                    poles_filter = [poles_filter; p_lower];
                    gain_filter = gain_filter*k_lower;
                end
                if(frequency_band_cutoffs(1) > 0) % 2nd order butterworth HP filter
                    [z_higher,p_higher,k_higher] = butter(2,frequency_band_cutoffs(1).*2./model.Probe.DAC.frequency, 'high');      % Doubled frequency, as "butter" requires cutt-off freq. normed by Nyquist rate
                    zeros_filter = [zeros_filter; z_higher];
                    poles_filter = [poles_filter; p_higher];
                    gain_filter = gain_filter*k_higher;
                end
                sos_filter = zp2sos(zeros_filter,poles_filter,gain_filter);
                filt = dfilt.df2sos(sos_filter);
                IR_filter_hat = freqz(filt, fft_length,'whole',model.Probe.DAC.frequency); % butterworth filter
                IR_filter_hat_zero_phase = abs(IR_filter_hat).^2; % zero-phase butterworth filter
            else
                fprintf('The map is computed for the following frequency band: [0,%d]\n',model.Probe.DAC.frequency/2);
                IR_filter_hat_zero_phase = ones(fft_length,1);
            end
            
            % Compute FFT of average or individual EIR
            if model.Probe.includeEir
                EIR_hat = fft(model.Probe.detector.EIR,fft_length,1)./fft_length.*IR_filter_hat_zero_phase;
            else
                EIR_hat = IR_filter_hat_zero_phase;
            end
            
            % Compute sensitivity for each detector
            for index_detector = number_detector
                
                for index_pixel = 1:model.Discretization.numberOfPixels
                    
                    % Compute FFT of TIR
                    if model.Probe.includeIndivEir
                        TPMR_sensitivity_hat = fft(model.Data.SPMR(:,model.Data.SPMRIndices(index_pixel,index_detector)),fft_length,1)./fft_length.*EIR_hat(:,ind_detector);
                    else
                        TPMR_sensitivity_hat = fft(model.Data.SPMR(:,model.Data.SPMRIndices(index_pixel,index_detector)),fft_length,1)./fft_length.*EIR_hat;
                    end
                    
                    % Integrate over TIR numerically
                    TPMR_sensitivity_map(index_pixel,index_detector) = trapz(abs(TPMR_sensitivity_hat).^2);
                end
            end
            
            % Average over all detectors and rescale to 1
            TPMR_sensitivity_map = mean(TPMR_sensitivity_map,2);
            TPMR_sensitivity_map = TPMR_sensitivity_map./max(abs(TPMR_sensitivity_map),[],'all');
            TPMR_sensitivity_map = reshape(TPMR_sensitivity_map,[model.Discretization.sizeOfPixelGrid(2),model.Discretization.sizeOfPixelGrid(1)]);
            
            % Plot map
            if length(number_detector) > 1
                figure_string = ['MSOT Probe ' model.Probe.probeId ' - Normalized Sensitivity Map For Detectors ' num2str(number_detector(1)) '_' num2str(number_detector(end)) '_freq_' num2str(frequency_band_cutoffs(1)/1e4) 'e4_' num2str(frequency_band_cutoffs(2)/1e4) 'e4'];
            else
                figure_string = ['MSOT Probe ' model.Probe.probeId ' - Normalized Sensitivity Map For Detector ' num2str(number_detector)  '_freq_' num2str(frequency_band_cutoffs(1)/1e4) 'e4_' num2str(frequency_band_cutoffs(2)/1e4) 'e4'];
            end
            figure('name',figure_string);
            set(gcf, 'Units', 'centimeters', 'Position', [0 1 15 12]);
            hold on;
            
            % Sensitivity map
            imagesc(model.Discretization.region.xPositionsFov(1,:),model.Discretization.region.zPositionsFov(:,1)',TPMR_sensitivity_map);
            
            % -3dB contour
            [~,c_TIR] = contour(model.Discretization.region.xPositionsFov,model.Discretization.region.zPositionsFov,TPMR_sensitivity_map,[0.5 0.5]);
            
            % Membrane
            if isa(model.Probe.membrane,'function_handle')
                plot(model.Discretization.region.xPositionsFov(1,:),model.Probe.membrane(model.Discretization.region.xPositionsFov(1,:)),'Color',[0, 0.4470, 0.7410]);
            end
            
            hold off;
            colormap gray;
            cb = colorbar;
            cb.TickLabelInterpreter = 'latex';
            caxis([0 1]);
            pbaspect([1 1 1]);
            xlim([model.Discretization.region.xPositionsFov(1,1) model.Discretization.region.xPositionsFov(1,end)]);
            ylim([model.Discretization.region.zPositionsFov(1,1) model.Discretization.region.zPositionsFov(end,1)]);
            c_TIR.LineColor = [0.8500, 0.3250, 0.0980];
            ylabel(cb, 'normalized TPMR energy','Interpreter','latex', 'FontSize',14);
            set(gca,'TickLabelInterpreter', 'latex', 'FontSize',14);
            xlabel('pixel position $x_{fov}$','Interpreter','Latex', 'FontSize',14);
            ylabel('pixel position $z_{fov}$','Interpreter','Latex', 'FontSize',14);
        end 
    end

    
    methods (Access=private)
        function sinograms = forward(model,initial_pressure)
            % Apply prefiltering
            coeffs_shift_invariant_space = model.Discretization.project_into_shift_invariant_space(initial_pressure);
            
            % Apply model matrix
            coeffs_shift_invariant_space = permute(reshape(coeffs_shift_invariant_space, model.Discretization.numberOfPixels, []), [2,1]);
            sinograms = forwardMex(coeffs_shift_invariant_space,...
                model.Data.SPMRIndices,...
                model.Data.absoluteTravelTimeIndices,...
                model.Data.SPMR,...
                model.Data.numTimeSamplesPerTransducer,...
                model.Data.SPMRTimestepIndexStartAndEnd,...
                model.Data.circularInterpolationRadiusInSamples);
            sinograms = permute(sinograms, [2,3,1]);
            
            % Convolve with impulse response
            sinograms = model.convolve_with_IR(sinograms);
        end
        
        function sinograms = forwardWithoutSIR(model,initial_pressure)
            % Apply prefiltering
            coeffs_shift_invariant_space = model.Discretization.project_into_shift_invariant_space(initial_pressure);
            
            % Apply model matrix
            coeffs_shift_invariant_space = permute(reshape(coeffs_shift_invariant_space, model.Discretization.numberOfPixels, []), [2,1]);
            sinograms = forwardDelayAndSumMex(coeffs_shift_invariant_space,...
                model.Data.absoluteTravelTimeIndices,...
                model.Data.numTimeSamplesPerTransducer,...
                model.Data.circularInterpolationRadiusInSamples);
            sinograms = permute(sinograms, [2,3,1]);
            
            % Convolve with impulse respones
            sinograms = model.convolve_with_IR(sinograms);
        end
        
        function transpose_coeffs_shift_invariant_space = transpose(model,sinograms)
            % Convolve with flipped impulse response
            sinograms = reshape(sinograms, model.Data.numTimeSamplesPerTransducer, model.Probe.detector.numOfTransducers, []);
            sinograms = model.transpose_convolve_with_IR(sinograms);
            
            % Apply transposed model matrix
            sinograms = permute(sinograms,[3,1,2]);        % permute matrix, so wavelength is first dimension (efficiency)
            transpose_pressure = transposeMex(sinograms,...
                model.Data.SPMRIndices,...
                model.Data.absoluteTravelTimeIndices,...
                model.Data.SPMR,...
                model.Data.SPMRTimestepIndexStartAndEnd,...
                model.Data.circularInterpolationRadiusInSamples);
            transpose_pressure = permute(transpose_pressure,[2,3,1]);               % permute matrix, so wavelength is last dimension (consistency)
            
            % Apply transpose prefiltering
            transpose_coeffs_shift_invariant_space = model.Discretization.transpose_project_into_shift_invariant_space(transpose_pressure);
        end
        
        function transpose_coeffs_shift_invariant_space = transposeWithoutSIR(model,sinograms)
            % Convolve with flipped impulse response
            sinograms = reshape(sinograms, model.Data.numTimeSamplesPerTransducer, model.Probe.detector.numOfTransducers, []);
            sinograms = model.transpose_convolve_with_IR(sinograms);
            
            % Apply transposed model matrix
            sinograms = permute(sinograms,[3,1,2]);       % permute matrix, so wavelength is first dimension (efficiency)
            transpose_pressure = transposeDelayAndSumMex(sinograms,...
                model.Data.absoluteTravelTimeIndices,...
                model.Data.circularInterpolationRadiusInSamples);
            transpose_pressure = permute(transpose_pressure,[2,3,1]);               % permute matrix, so wavelength is last dimension (consistency)
            
            % Apply transpose prefiltering
            transpose_coeffs_shift_invariant_space = model.Discretization.transpose_project_into_shift_invariant_space(transpose_pressure);
        end
        
        function sinograms = convolve_with_IR(model, sinograms)
            sinograms = ifft(fft(sinograms, size(model.Data.IRHat,1)) .* model.Data.IRHat);
            half_length_of_ir = floor(model.Probe.eirLength/2);
            sinograms = sinograms(half_length_of_ir+1 : end-half_length_of_ir, :, :);
        end
        
        function sinograms = transpose_convolve_with_IR(model, sinograms)
            sinograms = ifft(fft(sinograms, size(model.Data.IRFlippedHat,1)) .* model.Data.IRFlippedHat);
            half_length_of_ir = floor(model.Probe.eirLength/2);
            sinograms = sinograms(half_length_of_ir+1 : end-half_length_of_ir, :, :);
        end
        
        function funcs = generate_normalized_model_functions(model, forward, transpose, model_normalization_factor)
            if ~isempty(model_normalization_factor)
                fprintf('MSOTmodel: Use the normalization factor given as input to the model.\n');            
            else
                fprintf('MSOTmodel: Calculate greatest singular value of model to normalize it.\n');
                tic;
                unnormalizedModel = @(x,flag) combine_forward_and_transpose(x, flag, forward, transpose);
                model_normalization_factor = svds(unnormalizedModel, [model.Data.numTimeSamplesPerTransducer * model.Probe.detector.numOfTransducers, model.Discretization.numberOfPixels], 1, 'largest', 'MaxIterations', 50, 'Tolerance', 1e-8, 'FailureTreatment', 'keep');
                calcTime = toc;
                fprintf(['MSOTmodel: Greatest singular value calculated in ' num2str(calcTime) ' seconds.\n']);
            end
            
            fprintf(['MSOTmodel: The applied model normalization factor is ' num2str(model_normalization_factor) '\n']);
            funcs = [];
            funcs.applyForward = @(p0) forward(p0) / model_normalization_factor;
            funcs.applyTranspose = @(sigMat) transpose(sigMat) / model_normalization_factor;
        end
        
        function [hash, hashed_struct] = hash_to_load_and_save_acoustic_data(model)
            hashed_struct.coupling = model.Probe.coupling;
            hashed_struct.detector.x = model.Probe.detector.xPositionsOfTransducers;
            hashed_struct.detector.y = model.Probe.detector.zPositionsOfTransducers;
            hashed_struct.Medium = model.Medium;
            hashed_struct.region = model.Discretization.region;
            hashed_struct.useSingleSos = model.useSingleSpeedOfSound;
            
            % Replace membrane function handle with vector so that hashing is possible
            if model.useSingleSpeedOfSound
                hashed_struct.coupling.speedOfSound = hashed_struct.Medium.speedOfSound;
                hashed_struct.membrane = Inf*ones(size(model.Discretization.region.xPositionsFov(1,:)));
            else
                assert(isa(model.Probe.membrane,'function_handle'),'Dual speed of sound requires membrane model.');
                hashed_struct.membrane = model.Probe.membrane(model.Discretization.region.xPositionsFov(1,:));
            end
            
            % Hash data
            hash = DataHash(hashed_struct);
        end
        
        function [pixel_model_hash,pixel_model_string] = hash_to_load_and_save_SPMR(model)
            % Create string out of pixel model function handle for hashing
            pixel_model_string = func2str(model.Discretization.pixelModel);
            
            % Evaluate pixel model for hashing
            eval_bound = 10;
            pixel_model_eval = model.Discretization.pixelModel(-eval_bound:eval_bound);
            pixel_model_hash = DataHash(pixel_model_eval);
        end
        
        function [IR_hat,IR_flipped_hat] = get_impulse_response_for_model_with_spmr(model)
            IR = filter_butter_zero_phase(model.Probe.detector.EIR, model.Probe.DAC.frequency, [model.DataPreprocessing.filtCutoffMin, model.DataPreprocessing.filtCutoffMax], true);
            
            % Save IR in Fourier domain to enable fast convolution during model evaluations
            length_for_linear_fft = size(IR,1) + model.Data.numTimeSamplesPerTransducer-1;
            IR_hat = fft(IR, length_for_linear_fft);
            IR_flipped_hat = fft(flip(IR),  length_for_linear_fft);
        end
        
        function [IR_hat,IR_flipped_hat] = get_impulse_response_for_model_without_spmr(model)
            % select correct speed of sound
            if model.useSingleSpeedOfSound
                speed_of_sound_for_IR = model.Medium.speedOfSound;
            else
                speed_of_sound_for_IR = model.Probe.coupling.speedOfSound;
            end
            
            % convolve EIR and pixel model response
            time_samples_pixel_model_response = (-floor((model.Probe.eirLength-1)/2):floor((model.Probe.eirLength-1)/2))'/model.Probe.DAC.frequency;
            pixel_model_response = model.Discretization.pixelModel(time_samples_pixel_model_response*speed_of_sound_for_IR).*(-time_samples_pixel_model_response);  
            IR = convn(model.Probe.detector.EIR,pixel_model_response,'same');
            
            IR = filter_butter_zero_phase(IR, model.Probe.DAC.frequency, [model.DataPreprocessing.filtCutoffMin, model.DataPreprocessing.filtCutoffMax], true);
            
            % Save IR in Fourier domain to enable fast convolution during model evaluations
            length_for_linear_fft = size(IR,1) + model.Data.numTimeSamplesPerTransducer-1;
            IR_hat = fft(IR, length_for_linear_fft);
            IR_flipped_hat = fft(flip(IR),  length_for_linear_fft);
        end
        
        function [absolute_travel_time_indices,acoustic_data] = calculate_travel_time_indices(model)
            [hash, current_hashed_struct] = hash_to_load_and_save_acoustic_data(model);
            acoustic_file = which(['acoustic_file_' hash '.mat']);
            
            if not(isempty(acoustic_file))
                load(acoustic_file,'hashed_struct');
                assert(isequaln(hashed_struct, current_hashed_struct), 'Hash collision: Delete outdated travel time files');
                
                if model.includeSir
                    acoustic_data = load(acoustic_file,'travel_times','x_membrane_intersections');
                else
                    acoustic_data = load(acoustic_file,'travel_times');
                end
            else
                if model.includeSir || model.useSingleSpeedOfSound
                    [acoustic_data.x_membrane_intersections,acoustic_data.travel_times] = model.calculate_and_save_acoustic_data_for_model_with_SPMR_or_singSoS();
                else
                    acoustic_file_without_x_membrane_intersections = which(['acoustic_file_without_x_membrane_intersections_' hash '.mat']);
                    
                    if not(isempty(acoustic_file_without_x_membrane_intersections))
                        load(acoustic_file_without_x_membrane_intersections, 'hashed_struct');
                        assert(isequaln(hashed_struct, current_hashed_struct), 'Hash colision: Delete outdated travel time files');
                        acoustic_data = load(acoustic_file_without_x_membrane_intersections,'travel_times');
                        
                    else
                        acoustic_data.travel_times = calculate_and_save_acoustic_data_for_model_dualSoS_without_SPMR(model);
                    end
                end
            end
            
            indices_travel_time = acoustic_data.travel_times.*model.Probe.DAC.frequency - model.Probe.DAC.delayBeforeRecording;
            absolute_travel_time_indices = reshape(indices_travel_time, [], model.Probe.detector.numOfTransducers) - model.DataPreprocessing.numCroppedSamplesAtSinogramStart;
        end
        
        function travel_times = calculate_and_save_acoustic_data_for_model_dualSoS_without_SPMR(model)
            fprintf('Calculate travel times for dual speed of sound without SPMR correction.\n');
            resolution_x = (model.Discretization.region.xPositionsFov(end) - model.Discretization.region.xPositionsFov(1)) / (model.Discretization.sizeOfPixelGrid(1) - 1);
            resolution_z = (model.Discretization.region.zPositionsFov(end) - model.Discretization.region.zPositionsFov(1)) / (model.Discretization.sizeOfPixelGrid(2) - 1);
            
            assert(abs(resolution_x-resolution_z)<1e-12, 'Fast Marching Code non-applicable for non-uniform grid resolution. Please adjust Discretization.sizeOfPixelGrid or Discretization.region!');
            
            x_detector_on_grid_min = floor(min(model.Probe.detector.xPositionsOfTransducers)/resolution_x)*resolution_x;
            x_detector_on_grid_max = ceil(max(model.Probe.detector.xPositionsOfTransducers)/resolution_x)*resolution_x;
            z_detector_on_grid_min = floor(min(model.Probe.detector.zPositionsOfTransducers)/resolution_z)*resolution_z;
            z_detector_on_grid_max = ceil(max(model.Probe.detector.zPositionsOfTransducers)/resolution_z)*resolution_z;
            
            x_fov_min = min(min(model.Discretization.region.xPositionsFov(:)), x_detector_on_grid_min);
            x_fov_max = max(max(model.Discretization.region.xPositionsFov(:)), x_detector_on_grid_max);
            z_fov_min = min(min(model.Discretization.region.zPositionsFov(:)), z_detector_on_grid_min);
            z_fov_max = max(max(model.Discretization.region.zPositionsFov(:)), z_detector_on_grid_max);
            
            x_fov_vec = linspace(x_fov_min, x_fov_max, round((x_fov_max-x_fov_min)/resolution_x)+1);
            z_fov_vec = linspace(z_fov_min, z_fov_max, round((z_fov_max-z_fov_min)/resolution_z)+1);
            [x_positions_fov,z_positions_fov] = meshgrid(x_fov_vec,z_fov_vec);
            speed_of_sound_below_membrane = model.Medium.speedOfSound;
            speed_of_sound_above_membrane = model.Probe.coupling.speedOfSound;
            x_detectors = model.Probe.detector.xPositionsOfTransducers;
            z_detectors = model.Probe.detector.zPositionsOfTransducers;
            
            % Use fast marching code to calculate travel times
            tic;
            travel_times = calculate_travel_time_with_dual_speed_of_sound_without_intersec(model.Probe.membrane,x_positions_fov,z_positions_fov,x_detectors,z_detectors,speed_of_sound_below_membrane,speed_of_sound_above_membrane);
            toc;
            
            % Find indices in x/z that correspond to minimal/maximal entries of the x/z coordinates of the region Discretization
            [~, start_indices_x_of_fov_in_travel_time] = min(abs(x_positions_fov(1,:) - model.Discretization.region.xPositionsFov(1,1)));
            [~, end_indices_x_of_fov_in_travel_time] = min(abs(x_positions_fov(1,:) - model.Discretization.region.xPositionsFov(1,end)));
            
            [~, start_indices_z_of_fov_in_travel_time] = min(abs(z_positions_fov(:,1) - model.Discretization.region.zPositionsFov(1,1)));
            [~, end_indices_z_of_fov_in_travel_time] = min(abs(z_positions_fov(:,1) - model.Discretization.region.zPositionsFov(end,1)));
            
            travel_times = travel_times(start_indices_z_of_fov_in_travel_time:end_indices_z_of_fov_in_travel_time, start_indices_x_of_fov_in_travel_time:end_indices_x_of_fov_in_travel_time,:);
            
            assert(isequal(size(travel_times), [size(model.Discretization.region.xPositionsFov, 1) size(model.Discretization.region.xPositionsFov, 2) length(model.Probe.detector.xPositionsOfTransducers)]), 'Size of the travel time array does not match with Discretizationetization.');
            
            [hash, hashed_struct] = hash_to_load_and_save_acoustic_data(model);
            path_to_current_file = mfilename('fullpath');
            save(fullfile(path_to_current_file(1:end-10), '..', 'data', 'travel_time_files', ['acoustic_file_without_x_membrane_intersections_' hash '.mat']) , 'travel_times', 'hashed_struct'); 
        end
        
        function [x_membrane_intersections,travel_times] = calculate_and_save_acoustic_data_for_model_with_SPMR_or_singSoS(model)
            % Set model parameters for travel time computation
            if(model.useSingleSpeedOfSound)
                fprintf('MSOTmodel - Calculate travel times for single speed of sound.\n');
                membrane_model = @(x_mem) Inf;
                speed_of_sound_above_membrane = model.Medium.speedOfSound;
                speed_of_sound_below_membrane = model.Medium.speedOfSound;
            else
                fprintf('MSOTmodel - Calculate travel times for SPMR correction.\n');
                assert(isa(model.Probe.membrane,'function_handle'), 'Membrane has to be defined for dual speed of sound.');
                membrane_model = model.Probe.membrane;
                speed_of_sound_above_membrane = model.Probe.coupling.speedOfSound;
                speed_of_sound_below_membrane = model.Medium.speedOfSound;
            end

            % Compute travel time array and membrane intersections
            tic;
            [x_membrane_intersections,travel_times] = calculate_travel_time(membrane_model,...
                model.Probe.detector.radius,...
                model.Discretization.region.xPositionsFov,...
                model.Discretization.region.zPositionsFov,...
                model.Probe.detector.xPositionsOfTransducers,...
                model.Probe.detector.zPositionsOfTransducers,...
                speed_of_sound_below_membrane,...
                speed_of_sound_above_membrane);
            toc;
            
            % Create hash & store data
            [hash, hashed_struct] = hash_to_load_and_save_acoustic_data(model);
            path_to_current_file = mfilename('fullpath');
            save(fullfile(path_to_current_file(1:end-10), '..', 'data', 'travel_time_files', ['acoustic_file_' hash '.mat']), 'x_membrane_intersections', 'travel_times', 'hashed_struct');
        end
        
        function [SPMR,SPMR_timestep_index_start_and_end,SPMR_indices] = get_SPMR_and_SPMR_indices(model, acoustic_data)
            % Get polar coordinates of grid points in field of view
            [distances_from_transducers_to_virtual_fov, angles_azimuth_from_transducers_to_virtual_fov] = ...
                model.map_grid_point_coordinates_to_polar_coordinates_for_SPMR(acoustic_data);
            
            % Calculate SPMR or load precomputed
            [SPMR, SPMR_timestep_index_start_and_end, SPMR_min_distance, SPMR_stepsize_distance, SPMR_stepsize_angle_azimuth] = model.get_SPMR(distances_from_transducers_to_virtual_fov, angles_azimuth_from_transducers_to_virtual_fov);
            
            % Compute SPMR indices
            indices_in_SPMR_of_distances_from_transducers_to_virtual_fov = ...
                round(1 + (distances_from_transducers_to_virtual_fov-SPMR_min_distance) / SPMR_stepsize_distance); % '1+' due to matlab indexing
            indices_in_SPMR_angles_from_transducers_to_virtual_fov =...
                1 + round(angles_azimuth_from_transducers_to_virtual_fov*180/pi/SPMR_stepsize_angle_azimuth); % '1+' due to matlab indexing
            linearized_SPMR_indices = (indices_in_SPMR_of_distances_from_transducers_to_virtual_fov-1)*size(SPMR,2) ...
                + indices_in_SPMR_angles_from_transducers_to_virtual_fov;
            SPMR_indices =  int32(reshape(linearized_SPMR_indices, [], model.Probe.detector.numOfTransducers));
            
        end
        
        function [SPMR, SPMR_timestep_index_start_and_end, SPMR_min_distance, SPMR_stepsize_distance, SPMR_stepsize_angle_azimuth] = get_SPMR(model, distances_from_transducers_to_virtual_fov, angles_azimuth_from_transducers_to_virtual_fov)
            if model.useSingleSpeedOfSound
                speed_of_sound_for_SPMR = model.Medium.speedOfSound;
            else
                speed_of_sound_for_SPMR = model.Probe.coupling.speedOfSound;
            end

            [...
                required_min_distance_from_transducers_to_virtual_fov,...
                required_max_distance_from_transducers_to_virtual_fov,...
                required_distance_step_size_SPMR,...
                required_max_angle_azimuth_from_transducers_to_virtual_fov,...
                required_angle_azimuth_step_size_SPMR] = model.calculate_required_SPMR_distances_and_angels(distances_from_transducers_to_virtual_fov, angles_azimuth_from_transducers_to_virtual_fov,speed_of_sound_for_SPMR);
              
            % Check if precomputed polarSPMR is availabe and sufficient
            [pixel_model_hash,pixel_model_string] = model.hash_to_load_and_save_SPMR();
            SPMR_file = which(['polarSPMR_' model.Probe.probeId '_c_' num2str(speed_of_sound_for_SPMR) '_pixelModel_' num2str(pixel_model_hash) '.mat']);

            [...
                precomputed_SPMR_available,...
                precomputed_SPMR_version,...
                precomputed_SPMR_min_distance,...
                precomputed_SPMR_max_distance,...
                precomputed_SPMR_stepsize_distance,...
                precomputed_SPMR_max_angle_azimuth,...
                precomputed_SPMR_stepsize_angle_azimuth,...
                precomputed_SPMR_pixel_model_string] = model.load_properties_of_precomputed_polarSPMR(SPMR_file);
            
            sufficient_SPMR_file_available = ...
                precomputed_SPMR_available && ...
                strcmpi(precomputed_SPMR_version, model.VALID_SPMR_VERSION) && ...
                is_smaller_or_almost_equal(precomputed_SPMR_min_distance, required_min_distance_from_transducers_to_virtual_fov, model.SPMR_GRID_NUMERICAL_TOLERANCE) &&...
                is_larger_or_almost_equal(precomputed_SPMR_max_distance, required_max_distance_from_transducers_to_virtual_fov, model.SPMR_GRID_NUMERICAL_TOLERANCE) &&...
                is_smaller_or_almost_equal(precomputed_SPMR_stepsize_distance, required_distance_step_size_SPMR,  model.SPMR_GRID_NUMERICAL_TOLERANCE) &&...
                is_larger_or_almost_equal(precomputed_SPMR_max_angle_azimuth, required_max_angle_azimuth_from_transducers_to_virtual_fov,  model.SPMR_GRID_NUMERICAL_TOLERANCE) &&...
                is_smaller_or_almost_equal(precomputed_SPMR_stepsize_angle_azimuth, required_angle_azimuth_step_size_SPMR,   model.SPMR_GRID_NUMERICAL_TOLERANCE) &&...
                strcmpi(precomputed_SPMR_pixel_model_string,pixel_model_string);
            
                
            % Load or recompute SPMR
            if sufficient_SPMR_file_available
                fprintf('MSOTmodel: Precomputed SPMR for probe %s is available and sufficient for defined model. SPMR will be loaded.\n',model.Probe.probeId);
                load(SPMR_file, 'SPMR_polar','SPMR_timestep_index_start_and_end');               
                
                SPMR_min_distance = precomputed_SPMR_min_distance;
                SPMR_stepsize_distance = precomputed_SPMR_stepsize_distance;
                SPMR_stepsize_angle_azimuth = precomputed_SPMR_stepsize_angle_azimuth;
            else
                fprintf('MSOTmodel: Precomputed SPMR for probe %s is either not available or not sufficient for defined model. SPMR will be recomputed.\n',model.Probe.probeId);
                
                % If the precomputed SIR has the current version, select the more general parameters from "required_" and "precomputed_" to (re-)calculate the SIR 
                if strcmpi(precomputed_SPMR_version, model.VALID_SPMR_VERSION)
                    SPMR_min_distance = min(precomputed_SPMR_min_distance,required_min_distance_from_transducers_to_virtual_fov);
                    SPMR_max_distance = max(precomputed_SPMR_max_distance,required_max_distance_from_transducers_to_virtual_fov);
                    SPMR_stepsize_distance = min(precomputed_SPMR_stepsize_distance,required_distance_step_size_SPMR);
                    SPMR_max_angle_azimuth = max(precomputed_SPMR_max_angle_azimuth,required_max_angle_azimuth_from_transducers_to_virtual_fov);
                    SPMR_stepsize_angle_azimuth = min(precomputed_SPMR_stepsize_angle_azimuth,required_angle_azimuth_step_size_SPMR);
                else
                    SPMR_min_distance = required_min_distance_from_transducers_to_virtual_fov;
                    SPMR_max_distance = required_max_distance_from_transducers_to_virtual_fov;
                    SPMR_stepsize_distance = required_distance_step_size_SPMR;
                    SPMR_max_angle_azimuth = required_max_angle_azimuth_from_transducers_to_virtual_fov;
                    SPMR_stepsize_angle_azimuth  = required_angle_azimuth_step_size_SPMR;
                end
                
                [SPMR_polar,SPMR_timestep_index_start_and_end]  = calculate_polarSPMR(model.Probe.probeId,...
                	model.VALID_SPMR_VERSION,...
                    speed_of_sound_for_SPMR,...
                    model.Probe.detector.elevationRadiusOfSingleTransducer,...
                    model.Probe.detector.heightofSingleTransducer,...
                    model.Probe.detector.pitchOfTransducers,...
                    model.Probe.detector.seperationBetweenTransducers,...
                    model.Probe.DAC.frequency,...
                    SPMR_min_distance,...
                    SPMR_max_distance,...
                    SPMR_stepsize_distance,...
                    SPMR_max_angle_azimuth,...
                    SPMR_stepsize_angle_azimuth,...
                    0,...
                    1,...
                    model.Discretization.pixelModel,...
                    pixel_model_hash);
            end
            SPMR_polar = squeeze(SPMR_polar); 
            SPMR = permute(SPMR_polar,  [3 1 2]);
            SPMR_timestep_index_start_and_end = squeeze(SPMR_timestep_index_start_and_end);
            SPMR_timestep_index_start_and_end = permute(SPMR_timestep_index_start_and_end, [3 1 2]);
            SPMR_timestep_index_start_and_end = reshape(SPMR_timestep_index_start_and_end,2,size(SPMR,2)*size(SPMR,3));
            clear polarSPMR;
        end
        
        function [distances_from_transducers_to_virtual_fov,...
                angles_azimuth_from_transducers_to_virtual_fov] = map_grid_point_coordinates_to_polar_coordinates_for_SPMR(model,acoustic_data)
            % Transform grid point coordinates in field of view into polar coordinates for SPMR 
            x_transducers = reshape(model.Probe.detector.xPositionsOfTransducers,[1,1,model.Probe.detector.numOfTransducers]);
            z_transducers = reshape(model.Probe.detector.zPositionsOfTransducers,[1,1,model.Probe.detector.numOfTransducers]);
                
            if model.useSingleSpeedOfSound
                x_virtuals_fov_to_account_for_refraction = repmat(model.Discretization.region.xPositionsFov,[1,1,model.Probe.detector.numOfTransducers]);
                z_virtuals_fov_to_account_for_refraction = repmat(model.Discretization.region.zPositionsFov,[1,1,model.Probe.detector.numOfTransducers]);
            else
                speed_of_sound_for_SPMR = model.Probe.coupling.speedOfSound;
                x_membrane_intersections = acoustic_data.x_membrane_intersections;
                z_membrane_intersections = model.Probe.membrane(acoustic_data.x_membrane_intersections);
                x_virtuals_fov_to_account_for_refraction = x_transducers + (x_membrane_intersections - x_transducers)./sqrt((x_membrane_intersections - x_transducers).^2 + (z_membrane_intersections - z_transducers).^2).*acoustic_data.travel_times*speed_of_sound_for_SPMR;
                z_virtuals_fov_to_account_for_refraction = z_transducers + (z_membrane_intersections - z_transducers)./sqrt((x_membrane_intersections - x_transducers).^2 + (z_membrane_intersections - z_transducers).^2).*acoustic_data.travel_times*speed_of_sound_for_SPMR;
                nan_indices = isnan(x_membrane_intersections(:,:,1));
                
                % Note: If membrane intersections are not defined, set the virtual location to the actual location of the pixel
                x_virtuals_fov_to_account_for_refraction(isnan(x_membrane_intersections)) = repmat(model.Discretization.region.xPositionsFov(nan_indices),[length(x_transducers),1]);
                z_virtuals_fov_to_account_for_refraction(isnan(x_membrane_intersections)) = repmat(model.Discretization.region.zPositionsFov(nan_indices),[length(x_transducers),1]);
            end
            distances_from_transducers_to_virtual_fov = sqrt((x_virtuals_fov_to_account_for_refraction-x_transducers).^2 + (z_virtuals_fov_to_account_for_refraction-z_transducers).^2);
            angles_azimuth_from_transducers_to_virtual_fov = acos(((x_transducers-x_virtuals_fov_to_account_for_refraction).*x_transducers + (z_transducers-z_virtuals_fov_to_account_for_refraction).*z_transducers)./(sqrt(x_transducers.^2+z_transducers.^2).*distances_from_transducers_to_virtual_fov));
            
        end
        
        function [min_distance_from_transducers_to_virtual_fov,...
                max_distance_from_transducers_to_virtual_fov,...
                distance_step_size_SPMR_required,...
                max_angle_azimuth_from_transducers_to_virtual_fov,...
                angle_azimuth_step_size_SPMR_required] = calculate_required_SPMR_distances_and_angels(model,distances_from_transducers_to_virtual_fov,angles_azimuth_from_transducers_to_virtual_fov,speed_of_sound_for_SPMR)
            
            % Boundaries of SPMR grid in polar coordinates
            resolution_fov = min((model.Discretization.region.xPositionsFov(1,2)-model.Discretization.region.xPositionsFov(1,1)), (model.Discretization.region.zPositionsFov(2,1)-model.Discretization.region.zPositionsFov(1,1)));
            subsampling_factor_for_distance_step_size = 1450/speed_of_sound_for_SPMR; % Chosen so that for reasonable sos (1450 or higher), 'distance_step_size_SPMR_required' is at least resolution_fov.
            distance_step_size_SPMR_required = resolution_fov*subsampling_factor_for_distance_step_size; 
            
            min_distance_from_transducers_to_virtual_fov = min(distances_from_transducers_to_virtual_fov(:)) - distance_step_size_SPMR_required;
            max_distance_from_transducers_to_virtual_fov = max(distances_from_transducers_to_virtual_fov(:)) + distance_step_size_SPMR_required;
           
            angle_azimuth_step_size_SPMR_required = atan2(distance_step_size_SPMR_required, max_distance_from_transducers_to_virtual_fov)/pi*180;
            % Note: Min required angle is considered to be 0.
            max_angle_azimuth_from_transducers_to_virtual_fov = ceil((max(abs(angles_azimuth_from_transducers_to_virtual_fov(:)))/pi*180)*angle_azimuth_step_size_SPMR_required)/angle_azimuth_step_size_SPMR_required;
        end
        
        function [SPMR_file_available,...
                spmr_version,...
                min_distance_SPMR,...
                max_distance_SPMR,...
                stepsize_distance_SPMR,...
                max_angle_azimuth_SPMR,...
                stepsize_angle_azimuth_SPMR,...
                pixel_model_string_SPMR] = load_properties_of_precomputed_polarSPMR(~,SPMR_file)
            SPMR_file_available = ~isempty(SPMR_file);
            
            if SPMR_file_available
                % Load SPMR version
                mat_file_var = who('-file',SPMR_file,'spmr_version');
                if ~isempty(mat_file_var)
                    load(SPMR_file,mat_file_var{1});
                else
                    spmr_version = 0;
                end
                
                % Load distance and azimuth angle matrices of DIR
                load(SPMR_file, 'angle_azimuth_mat_SPMR', 'distance_mat_SPMR','pixel_model_string_SPMR');
                min_distance_SPMR = distance_mat_SPMR(1);
                max_distance_SPMR = distance_mat_SPMR(end);
                stepsize_distance_SPMR = distance_mat_SPMR(2)-distance_mat_SPMR(1);
                max_angle_azimuth_SPMR = angle_azimuth_mat_SPMR(end);
                stepsize_angle_azimuth_SPMR = angle_azimuth_mat_SPMR(2)-angle_azimuth_mat_SPMR(1);
            else
                spmr_version = 0;
                min_distance_SPMR = NaN;
                max_distance_SPMR = NaN;
                stepsize_distance_SPMR = NaN;
                max_angle_azimuth_SPMR = NaN;
                stepsize_angle_azimuth_SPMR = NaN;
                pixel_model_string_SPMR = '';
            end
        end
    end
end

