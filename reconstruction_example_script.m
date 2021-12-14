run('startup_reconstruction.m');

%% Define model
device_probe_id = 'example_system';               % see Probe.m for all available probes
use_eir = true;
use_indiv_eir = true;
use_sir = true;
use_single_speed_of_sound = true;
num_cropped_samples_at_sinogram_start = 1;
filt_cutoff_min = 1e5;
filt_cutoff_max = 12e6;
field_of_view = [-0.02 0.02 -0.02 0.02];          % [x_fov_min x_fov_max z_fov_min z_fov_max]
number_of_grid_points_fov = [401 401];            % [grid_points_x_dimension grid_points_z_dimension]
speed_of_sound_tissue = 1465;
model_normalization_factor = [];                  % if empty, the model is normalized so that its largest singular values is 1.

model = define_model_for_reconstruction(field_of_view, number_of_grid_points_fov, device_probe_id, use_eir, use_indiv_eir, use_sir, use_single_speed_of_sound, speed_of_sound_tissue, num_cropped_samples_at_sinogram_start, filt_cutoff_min, filt_cutoff_max, model_normalization_factor);

%% Execute reconstruction
% Load  data
data_raw = load_raw_data_for_reconstruction(model); % TODO: Replace the placeholder implementation in this function to load data from real measurements

% Preprocess the signals
data_precrop = crop_first_n_signals(data_raw,  model.DataPreprocessing.numCroppedSamplesAtSinogramStart);
data_precrop_windowed = apply_butterworth_window_to_sinogram(data_precrop, 2, 300, size(data_precrop,1)-200);
data_filt = filter_butter_zero_phase(data_precrop_windowed, model.Probe.DAC.frequency, [model.DataPreprocessing.filtCutoffMin, model.DataPreprocessing.filtCutoffMax],true);
broken_transducers = [];
data_filt = interpolate_signals_of_broken_transducers(data_filt, broken_transducers);

% Reconstruction and save results
save_dir = 'path/to/output/directory';
file_name = 'example_reconstruction';

% Non-negative model-based reconstruction with Shearlet L1 regularization
lambda = 1e-4;
num_iterations = 50;
rec_img_nn = rec_nn_with_Shearlet_reg(model, data_filt, num_iterations, lambda);

% Save results
niftiwrite(rec_img_nn, fullfile(save_dir, [file_name '.nii']));


