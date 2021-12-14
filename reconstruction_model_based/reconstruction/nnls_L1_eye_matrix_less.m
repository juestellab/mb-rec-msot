function [rec_img, objective_per_iter, execution_time_per_iter] = nnls_L1_eye_matrix_less(model_func, target_vector, initial_img, etaReg, lambdaParam, nIter, cgIter, calculate_objective_per_iter)
% performing NNLS with L1 regularization and identity reg. matrix
% uses method from: "A Split Bregman Method for non-negative Sparsity
%                   Penalized Least Sq."
% Inputs: 
% > model_func:    function for computing the forward model
% > target_vector: vector of transducer data (we want the forward of rec_img to match them)
% > initial_img:   initial image - guess of the solution (usually set as vector 0)
%                  the x0 of the optimization method
% > etaReg:        eta regularization param in the method
% > lambdaRel:     lambda param in the method (usually small)
% > nIter:         max. number of outer iteration
% > cgIter:        max. number of inner iteration (cgs iterations)
% > calculate_objective_per_iter:   Flag if objective should be
%                                   calculated at the end of every outer
%                                   iteration (one additional model evaluation)
%
% Outputs: 
% > rec_img:                 solution - reconstructed image
% > objective_per_iter:      vector containing the objective value f(x) = ||Ax-b||^2 + eta*|u| 
%                            after each outer iteration
% > execution_time_per_iter: vector containing the execution time after
%                            each outer iteration (in sec)
% terminates if nIter is reached or epsStepRel condition is met

% initialization (to zero)
b_param = zeros(size(initial_img)); % b parameter of method
img_inner_iter = initial_img; % u param of method
% initialize f values and execution times
objective_per_iter = zeros(nIter, 1);
execution_time_per_iter = zeros(nIter, 1);
% define function handler for computing (lambda * A' * A + I)*u
L1_method_operator = @(v) lambdaParam * model_func( model_func(v, 'notransp'), 'transp' ) + v;
inner_tol = 1e-8;

% outter iterations
tic;
for i = 1:nIter
    start_time = toc;
    % update d parameter
    d_param = P(img_inner_iter - b_param);
    % update img_inner_iter: solve inner system of equations to get update
    b_inner = P(d_param) + b_param + lambdaParam * model_func(target_vector, 'transp') - lambdaParam * etaReg * ones(size(b_param));
    % solve system with cgs and get uu update (inner iteration)
    [img_inner_iter, ~] = cgs(L1_method_operator, b_inner, inner_tol, cgIter);
    % update b_param
    b_param = b_param + P(d_param) - img_inner_iter;
    
    % iteration end, store values
    if calculate_objective_per_iter 
        objective_per_iter(i) = norm( model_func(img_inner_iter, 'notransp') - target_vector )^2 + etaReg * sum (abs(img_inner_iter));
    end
    execution_time_per_iter(i) = toc - start_time;
end

% done, return solution
rec_img = img_inner_iter;
end

function u = P(d)
% projection operation: for input vector d, returns d(i) if d(i) > 0,
% or replaces d(i) with 0 otherwise
u = (d > 0) .* d;

end


