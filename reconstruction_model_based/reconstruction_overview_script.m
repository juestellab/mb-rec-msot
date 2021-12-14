run('../startup_reconstruction.m');

%% Define model
device_probe_id = 'example_system';                   % see Probe.m for all available probes
use_eir = true;
use_indiv_eir = true;
use_sir = false;
use_single_speed_of_sound = true;
num_cropped_samples_at_sinogram_start = 1;
filt_cutoff_min = 1e5;
filt_cutoff_max = 12e6;
field_of_view = [-0.02 0.02 -0.02 0.02];              % [x_fov_min x_fov_max y_fov_min y_fov_max]
number_of_grid_points_fov = [101 101];                % [grid_points_x_dimension grid_points_y_dimension]
speed_of_sound_tissue = 1465;

model = define_model_for_reconstruction(field_of_view, number_of_grid_points_fov, device_probe_id, use_eir, use_indiv_eir, use_sir, use_single_speed_of_sound, speed_of_sound_tissue, num_cropped_samples_at_sinogram_start, filt_cutoff_min, filt_cutoff_max);

%% Create test image (initial pressure distribution)
test_img = 0.01*ones(model.Discretization.sizeOfPixelGrid);
test_img(model.Discretization.region.zPositionsFov<0.015) = 0.1;
r = sqrt(model.Discretization.region.xPositionsFov.^2+model.Discretization.region.zPositionsFov.^2);
test_img(r<0.001) = 1;
test_img(10:20,23:30) = 1;

% Comment-in next line to simulate a frame consisting of 28 wavelengths
%test_img = repmat(test_img,1,1,28);

%% Simulate recorded pressure signals for the defined test image with the forward model
tic;
test_sig = model.Funcs.applyForward(test_img);
t_single_forward = toc;
test_sig = test_sig + 0.01 * randn(size(test_sig));
test_sig = reshape(test_sig,model.Probe.detector.numOfTransducers,size(test_img,3));

%% Obtain backprojection-like reconstruction with the transpose model
tic;
transp_img = model.Funcs.applyTranspose(test_sig);
t_single_transpose = toc;
transp_img = reshape(transp_img, model.Discretization.sizeOfPixelGrid(1), model.Discretization.sizeOfPixelGrid(2), []);

%% Model-based reconstructions
% Non-negative reconstruction without regularization
[rec_img_nn2, lCurveErrImg, lCurveErrReg, lCurveErrReg_2, err, t] = rec_nn_with_L2_reg(model, test_sig, 400, 0, [], [], 0, [], []);
figure; imagesc(rec_img_nn2)

% Non-negative reconstruction with L2 and L2 Laplace regularization
RegL2 = @(x) x;
RegL2T = @(x) x;
lambdaL2 = 5e-4;
RegL2_2 = @(x) laplacian_per_wavelength(reshape(x, model.Discretization.sizeOfPixelGrid(2), model.Discretization.sizeOfPixelGrid(1), []));
RegL2T_2 = @(x) laplacian_per_wavelength(reshape(x, model.Discretization.sizeOfPixelGrid(2), model.Discretization.sizeOfPixelGrid(1), []));
lambdaL2_2 = 5e-3;
lambdaL1 = 1e-5;
[rec_img_nnReg2, lCurveErrImg, lCurveErrReg, lCurveErrReg_2, err, t] = rec_nn_with_L2_reg(model, test_sig, 50, lambdaL2, RegL2, RegL2T, lambdaL2_2, RegL2_2, RegL2T_2); %#ok<*ASGLU>
figure; imagesc(rec_img_nnReg2)

% Non-negative reconstruction with L1 eye reg. matrix
[rec_img_nnReg_L1,  lCurveErrImg, lCurveErrReg, err, tReg] = rec_nn_with_L1_eye_reg(model, test_sig, 50, 1e-3);
figure; imagesc(rec_img_nnReg_L1)

% TV non-negative limited view reconstruction
[rec_img_nnReg_TV, lCurveErrImg, lCurveErrReg, err, t] = rec_nn_with_TV_reg(model, test_sig, 50, 1e-3);
figure; imagesc(rec_img_nnReg_TV)

% Shearlet non-negative limited view reconstruction
[rec_img_nnReg_Shearlet, lCurveErrImg, lCurveErrReg, err, t] = rec_nn_with_Shearlet_reg(model, test_sig, 50, 1e-4);
figure; imagesc(rec_img_nnReg_Shearlet)

