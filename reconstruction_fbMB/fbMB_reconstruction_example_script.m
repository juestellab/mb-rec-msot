run(fullfile('..', 'startup_reconstruction.m'));

%% Define model
device_probe_id = 'example_system';
use_eir = true;
use_indiv_eir = true;
use_sir = false;
use_single_speed_of_sound = true;
num_cropped_samples_at_sinogram_start = 1;
filt_cutoff_min = 1e5;
filt_cutoff_max = 12e6;
field_of_view = [-0.02 0.02 -0.02 0.02];
number_of_grid_points_fov = [401 401];
speed_of_sound_tissue = 1465;
model_normalization_factor = [];

model = define_model_for_reconstruction(field_of_view, number_of_grid_points_fov, device_probe_id, use_eir, use_indiv_eir, use_sir, use_single_speed_of_sound, speed_of_sound_tissue, num_cropped_samples_at_sinogram_start, filt_cutoff_min, filt_cutoff_max, model_normalization_factor);

%% Define frequency decomposition model
% filters
F1_min = 1e5;
F1_max = 2e6;
F{1} = @(sigMat) filter_butter_zero_phase(sigMat,model.Probe.DAC.frequency,[F1_min F1_max],0);

F2_min = 1e5;
F2_max = 13e6;
F{2} = @(sigMat) filter_butter_zero_phase(sigMat,model.Probe.DAC.frequency,[F2_min F2_max],0);

F3_min = 7e6;
F3_max = 13e6;
F{3} = @(sigMat) filter_butter_zero_phase(sigMat,model.Probe.DAC.frequency,[F3_min F3_max],0);

% lambdas
lambdaF = 1*[1,1,1];
lambdaR = 0*[1,1,1];
lambda0 = 2e-1;

% regularization
R0 = @(x) x;
R0T = @(x) x;
R = {[],[],[]};
RT = {[],[],[]};

[s,Ns] = freq_decomp_model(vec(zeros(number_of_grid_points_fov(1),number_of_grid_points_fov(2),length(F))),model.Funcs.applyForward,R0,F,R,lambda0,lambdaF,lambdaR);

% Define fbMB model
fbMB_forward = @(x) freq_decomp_model(x,model.Funcs.applyForward,R0,F,R,lambda0,lambdaF,lambdaR);
fbMB_transpose = @(s) freq_decomp_transp_model(s,model.Funcs.applyTranspose,R0T,F,RT,lambda0,lambdaF,lambdaR,Ns);

model_fbMB = @(x,flag) combine_forward_and_transpose(x, flag, fbMB_forward, fbMB_transpose);

%% Test reconstruction using a numerical phantom
% Define numerical phatnom
vx = linspace(field_of_view(1),field_of_view(2),number_of_grid_points_fov(1));
vy = linspace(field_of_view(1),field_of_view(2),number_of_grid_points_fov(2));
[x,y] = meshgrid(vx,vy);
r = sqrt(x.^2+y.^2);

test_img = zeros(size(r));
for j = 1:10
    dr = rand(1)*0.003;
    Rr = rand(1)*0.008;
    xr = 0.0125*rand(1) - 0.00625;
    yr = 0.0125*rand(1) - 0.00625;
    rr = sqrt((x-xr).^2+(y-yr).^2);
    test_img = test_img + rand(1)*(rr <= Rr).*(rr >= Rr-dr);
end
test_sig = model.Funcs.applyForward(test_img);
%test_sig = awgn(test_sig, 30); % Communications toolbox required

% Run fbMB reconstruction
N_inner = 5;
N_outer = 40;
target_vector = [vec(test_sig); zeros(size(test_sig(:))); zeros(size(test_sig(:))); zeros(size(test_sig(:))); zeros(size(test_img(:)))];
initial_guess = model_fbMB(target_vector,'transp');
lambda = 1;
beta = 1e-3;
calculate_objective_per_iter = 0;
[Coeffs, err, t] = nnls_L2_less(model_fbMB, target_vector, initial_guess, lambda, beta, N_outer, N_inner, calculate_objective_per_iter);

N = 401;
img = zeros(N,N,length(F));
for j = 1:length(F)
    img(:,:,j) = rot90(reshape(Coeffs(((j-1)*N^2+1):(j*N^2)),N,[]),2);
end
img_freq = sum(img,3);

% Run reference model-based reconstruction
lambdaL2 = 5e-4;
[img_ref, lCurveErrImg, lCurveErrReg, Coeffs_ref, errReg, tReg] = rec_nn_with_L2_reg(model, test_sig, N_outer,lambdaL2,R0, R0T);

% Visualization
Mfull = max(vec([test_img,img_freq]));
Mfreq = max(vec(img));

figure(1);
subplot(2,3,1); imagesc(vx,vy,rot90(test_img,2)); axis equal tight; colorbar; title('ground truth'); caxis([0 Mfull]);
subplot(2,3,2); imagesc(vx,vy,img_ref); axis equal tight; colorbar; title('reference reconstruction'); caxis([0 Mfull]);
subplot(2,3,3); imagesc(vx,vy,img_freq); axis equal tight; colorbar; title('freq. decomp. reconstruction'); caxis([0 Mfull]);
subplot(2,3,4); imagesc(vx,vy,img(:,:,1)); axis equal tight; colorbar; title([num2str(F1_min/1e6) ' to ' num2str(F1_max/1e6) ' MHz part']); caxis([0 Mfreq]);
subplot(2,3,5); imagesc(vx,vy,img(:,:,2)); axis equal tight; colorbar; title([num2str(F2_min/1e6) ' to ' num2str(F2_max/1e6) ' MHz part']); caxis([0 Mfreq]);
subplot(2,3,6); imagesc(vx,vy,img(:,:,3)); axis equal tight; colorbar; title([num2str(F3_min/1e6) ' to ' num2str(F3_max/1e6) ' MHz part']); caxis([0 Mfreq]);

sig1 = model.Funcs.applyForward(img(:,:,1));
sig2 = model.Funcs.applyForward(img(:,:,2));
sig3 = model.Funcs.applyForward(img(:,:,3));
sig_freq = model.Funcs.applyForward(img_freq);
sig_ref = model.Funcs.applyForward(img_ref);

Sfull = max(vec(log(abs([fft(model.Funcs.applyForward(test_img)),fft(sig_freq),fft(sig_ref)]))));
Sfreq = max(vec(log(abs([fft(sig1),fft(sig2),fft(sig3)]))));

L = size(test_sig,1);
Fs = model.Probe.DAC.frequency;
fs = Fs*linspace(-L/2,L/2,L)/L/1e6; % MHz
ds = 1:size(test_sig,2);

figure(2);
subplot(2,3,1); imagesc(ds,fs,fftshift(log(abs(fft(model.Funcs.applyForward(test_img)))))); axis tight; colorbar; title('ground truth'); caxis([-10 Sfull]); ylim([0 max(fs)]);
subplot(2,3,2); imagesc(ds,fs,fftshift(log(abs(fft(sig_ref))))); axis tight; colorbar; title('reference reconstruction'); caxis([-10 Sfull]); ylim([0 max(fs)]);
subplot(2,3,3); imagesc(ds,fs,fftshift(log(abs(fft(sig_freq))))); axis tight; colorbar; title('freq. decomp. reconstruction'); caxis([-10 Sfull]); ylim([0 max(fs)]);
subplot(2,3,4); imagesc(ds,fs,fftshift(log(abs(fft(sig1))))); axis tight; colorbar; title([num2str(F1_min/1e6) ' to ' num2str(F1_max/1e6) ' MHz part']); caxis([-10 Sfreq]); ylim([0 max(fs)]);
subplot(2,3,5); imagesc(ds,fs,fftshift(log(abs(fft(sig2))))); axis tight; colorbar; title([num2str(F2_min/1e6) ' to ' num2str(F2_max/1e6) ' MHz part']); caxis([-10 Sfreq]); ylim([0 max(fs)]);
subplot(2,3,6); imagesc(ds,fs,fftshift(log(abs(fft(sig3))))); axis tight; colorbar; title([num2str(F3_min/1e6) ' to ' num2str(F3_max/1e6) ' MHz part']); caxis([-10 Sfreq]); ylim([0 max(fs)]);
