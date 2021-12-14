function [rec_img, lCurveErrImg, lCurveErrReg, err, t, shearletTrafo] = rec_nn_with_Shearlet_reg(model, data_filt, num_iter, regParam)
% Reconstruction method for Shearlet regularization
% Adaptation from paper: "Total-Variation Regularization With Bound Constraints"
% wrapper method for function: "nnls_sparsa_based_reg_matrix_less" 
% Inputs:
% > model:     the optoacoustic model
% > data_filt: the transducer data (num_sensors x num_samples x num_wavelengths)
% > num_iter:  number of iterations of optimization method
% > regParam:  strength of applied regularization. Suggested value: 1.6e-3
% Outputs:
% > rec_img:      the reconstructed image (high x width) or a stack of images
%                 (high x width x num_wavelengths), in case of multiple wavelengths
% > lCurveErrImg: the reconstruction error (L2 norm between data_filt and
%   the samples   that rec_img would give us
% > lCurveErrReg: regularization value on reconstructed image
% > err:          vector of objective values (reconstruction + regularization
%                 errors) for each iteration
% > t:            execution time for each iteration (in sec)


% get model size 
N = model.Discretization.sizeOfPixelGrid; % img width and height [Nx, Nz]

% set up shearlet system (see related methods)
shearletSystem = SLgetShearletSystem2D(0, N(2), N(1), 4);
shearletTrafo = @(x) SLsheardec2D(reshape(x, N(2), N(1)), shearletSystem);
shearletRecon = @(c) vec(SLshearrec2D(c, shearletSystem));

% set params to recommended values
num_iter_nn_inner = 2;
num_iter_nn_outer = num_iter; % typically 200 - 300
lambda = 1;
beta = 1e-3;
calculate_objective_per_iter = false;
Psi = @(x,tau) PsiForShearletReg(reshape(x, N(2), N(1), []), tau, shearletTrafo, shearletRecon);
Phi = @(x) PhiForShearletReg(reshape(x, N(2), N(1), []), shearletTrafo);

% build model function and target vector
applyModel = @(x,flag) combine_forward_and_transpose(x, flag, model.Funcs.applyForward, model.Funcs.applyTranspose);
targetVector = data_filt(:);

% set transpose as start value
transpose_as_start_value = applyModel(targetVector, 'transp');

% run Shearlet reconstruction method 
tic;
[rec_img, err, t] = nnls_sparsa_based_reg_matrix_less(applyModel, Psi, Phi, targetVector, transpose_as_start_value, lambda, beta, regParam, num_iter_nn_outer, num_iter_nn_inner, calculate_objective_per_iter);

% execution time
calcTime = toc;
fprintf(['Finished nn reconstruction in ' num2str(calcTime) ' seconds.\n']);
% data error
lCurveErrImg = norm(vec(model.Funcs.applyForward(rec_img) - data_filt),2)^2;
% regularization error
lCurveErrReg = Phi(rec_img);
% process reconstructed image (rotate, etc.)
rec_img = reshape(rec_img, N(2), N(1), []);

end
