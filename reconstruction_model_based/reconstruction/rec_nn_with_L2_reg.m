function [rec_img, lCurveErrImg, lCurveErrReg, lCurveErrReg_2, err, t] = rec_nn_with_L2_reg(model, data_filt, num_iter_nn_outer, lambdaL2, RegL2, RegL2T, lambdaL2_2, RegL2_2, RegL2T_2)
% Reconstruction method for L2 non-negative regularization
% Adaptation from paper: "Total-Variation Regularization With Bound Constraints"
% wrapper method for function: "nnls_L2_less" 
% Inputs:
% > model:      the optoacoustic model
% > data_filt:  the transducer data (num_sensors x num_samples x num_wavelengths)
% > num_iter:   number of iterations of optimization method
% > lambdaL2:   strength of L2 regularization. Suggested value: 5e-4
% > RegL2:      function computing the L2 regularization term (see rec. overview) 
% > RegL2T:     function computing the transpose of RegL2
% > lambdaL2_2: strength of Laplace regularization. Suggested value: 5e-3
% > RegL2:      function computing the Laplace regularization term (see rec. overview) 
% > RegL2T:     function computing the transpose of RegL2T
% Outputs:
% > rec_img:      the reconstructed image (high x width) or a stack of images
%                 (high x width x num_wavelengths), in case of multiple wavelengths
% > lCurveErrImg: the reconstruction error (L2 norm between data_filt and
%   the samples   that rec_img would give us
% > lCurveErrReg: regularization value on reconstructed image
% > err:          vector of objective values (reconstruction + regularization
%                 errors) for each iteration
% > t:            execution time for each iteration (in sec)
if nargin < 9
    lambdaL2_2 = [];
    RegL2_2 = [];
    RegL2T_2 = [];
end

% get model size - set parameters to recommended values
N = model.Discretization.sizeOfPixelGrid; % img width and height [Nx, Nz]
num_iter_nn_inner = 2;
lambda = 1;
beta = 1e-3;
calculate_objective_per_iter = false;

% build model function and target vector
% incoorporate reg. functions into model
[applyModel, targetVector] = incorporate_l2_reg_into_model_functions(data_filt, N, model.Funcs.applyForward, model.Funcs.applyTranspose, lambdaL2, RegL2, RegL2T, lambdaL2_2, RegL2_2, RegL2T_2);
% set transpose as start value
transpose_as_start_value = applyModel(targetVector, 'transp');

% run reconstruction
tic;
[rec_img, err, t] = nnls_L2_less(applyModel, targetVector, transpose_as_start_value, lambda, beta, num_iter_nn_outer, num_iter_nn_inner, calculate_objective_per_iter);

% execution time
calcTime = toc;
fprintf(['Finished nn reconstruction in ' num2str(calcTime) ' seconds.\n']);
% calculate values and errors
lCurveErrImg = norm(vec(model.Funcs.applyForward(rec_img) - data_filt),2)^2;
% L2 reg error
if isa(RegL2,'function_handle')
    lCurveErrReg = norm(vec(RegL2(rec_img)),2).^2;
else
    lCurveErrReg = 0.0;
end
% L2_2 reg error
if isa(RegL2_2,'function_handle')
    lCurveErrReg_2 = norm(vec(RegL2_2(rec_img)),2).^2;
else
    lCurveErrReg_2 = 0.0;
end
% process reconstructed image (rotate, etc.)
rec_img = reshape(rec_img, N(2), N(1), []);
end
