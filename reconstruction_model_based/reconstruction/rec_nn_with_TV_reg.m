function [rec_img, lCurveErrImg, lCurveErrReg, err, t] = rec_nn_with_TV_reg(model, data_filt, num_iter, regParam)
% Reconstruction method for TV regularization
% from paper: "Total-Variation Regularization With Bound Constraints"
% wrapper method for function "nnls_sparsa_based_reg_matrix_less" 
% Inputs:
% > model:     the optoacoustic model
% > data_filt: the transducer data (num_sensors x num_samples x num_wavelengths)
% > num_iter:  number of iterations of optimization method
% > regParam:  strength of applied regularization. Suggested value: 2.5e-2
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
N = model.Discretization.sizeOfPixelGrid;
% note: this TV implementation required squared images because of 'DiffOper'
if N(1) ~= N(2)
    error('This TV implementation works only with square images!');
end
%N = N(1);
num_iter_nn_inner = 2;
num_iter_nn_outer = num_iter;
lambda = 1;
beta = 1e-3;
calculate_objective_per_iter = false;

% build model function and target vector
applyModel = @(x,flag) combine_forward_and_transpose(x, flag, model.Funcs.applyForward, model.Funcs.applyTranspose);
targetVector = data_filt(:);

% set transpose as start value
transpose_as_start_value = applyModel(targetVector, 'transp');

% define Psi and Phi functions for TV
[B, Bt, BtB] = DiffOper(N(1)); % DiffOper requires squared images.
Psi = @(x,tau) PsiForTVReg(reshape(x,N(2),N(1),[]),tau,B,Bt,BtB);
Phi = @(x) PhiForTVReg(reshape(x, N(2), N(1), [])); % TV norm term

% run TV reconstruction method 
tic;
[rec_img, err, t] = nnls_sparsa_based_reg_matrix_less(applyModel, Psi, Phi, targetVector, transpose_as_start_value, lambda, beta, regParam, num_iter_nn_outer, num_iter_nn_inner, calculate_objective_per_iter);

calcTime = toc;
fprintf(['Finished nn reconstruction in ' num2str(calcTime) ' seconds.\n']);
lCurveErrImg = norm(vec(model.Funcs.applyForward(rec_img) - data_filt),2)^2;
lCurveErrReg = Phi(rec_img);
rec_img = reshape(rec_img, N(2), N(1), []);

end
