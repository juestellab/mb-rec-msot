function [rec_img, objective_per_iter, execution_time_per_iter] = nnls_L2_less(model_func, target_vector, initial_img, lambda, beta, nIter, nIter_inner, calculate_objective_per_iter)
% performing NNLS with L2 regularization
% adapts method from: "Total-Variation Regularization With Bound Constraints"
% Note: this method was originally derived by the authors for the TV norm,
% but we adapt it and use it also for L2.
% Inputs: 
% > model_func:    function for computing the forward model
% > target_vector: vector of transducer data (we want the forward of rec_img to match them)
% > initial_img:   initial image - guess of the solution (usually set as vector 0)
%                  the x0 of the optimization method
% > lambda:        lambda param of the method
% > beta:          beta param of the method
% > nIter:         max. number of outer iteration
% > nIter_inner:   max. number of inner iteration (usually set below 5, 
%                  depending on the problem size)
% > calculate_objective_per_iter:   Flag if objective should be
%                                   calculated at the end of every outer
%                                   iteration (one additional model evaluation)
% Outputs: 
% > rec_img:                 solution - reconstructed image
% > objective_per_iter:      vector containing the objective value f(x) = ||Ax-b||^2 + reg(u) 
%                            after each outer iteration
% > execution_time_per_iter: vector containing the execution time after
%                            each outer iteration (in sec)
% terminates if nIter is reached or epsStepRel condition is met

% initialization
img_out_iter = zeros(size(initial_img));
lagrange_multip_for_convergence = img_out_iter;
num_pixels = length(initial_img);
num_target_values = length(target_vector);
inner_tol = 1e-8; % inner iteration tolerance for SPARSA

% initialize objective values and execution times
objective_per_iter = zeros(nIter, 1);
execution_time_per_iter = zeros(nIter, 1);

% define function to apply model and w target (augmented matrix)
% according to the method, we need to solve the following problem in each
% inner iteration:
% (lambda/2) * ||Au - b||^2 + (beta/2) * ||w - u||^2
% We can write this as:
% 1/2 * ||A_aug * u - b_aug||^2,
% where: A_aug = [sqrt(lambda) * A; sqrt(beta) * I] (augmented model matrix) 
%        b_aug = [sqrt(lambda) * b; sqrt(beta) * w] (augm. target vector)
% this function is needed for lsqr
linear_operator_augmented = @(inp, flag) add_nn_constraint_term_to_model_matrix(inp, flag, model_func, lambda, beta, num_pixels, num_target_values);

% outter iterations
tic;
for i = 1:nIter
    start_time = toc;
    
    % step 1: solve augmented problem with LSQR
    % apply LSQR to get u update (nIter_inner)
    [img_inner_iter,~,~,~,~,~] = lsqr(linear_operator_augmented, [sqrt(lambda)*target_vector', sqrt(beta)*(img_out_iter+lagrange_multip_for_convergence)']', inner_tol, nIter_inner, [], [], img_out_iter);
    
    % step 2: solve 10 from paper (w update, this is simply clipping)
    % this is the image - solution of the outer iteration
    img_out_iter = max(img_inner_iter - lagrange_multip_for_convergence, 0);
    
    % step 3: v update (see paper)
    % this is a Lagrangian multiplier ensuring convergence
    lagrange_multip_for_convergence = lagrange_multip_for_convergence + img_out_iter - img_inner_iter;
    
    % iteration end, store values
    if calculate_objective_per_iter 
        objective_per_iter(i) = (lambda / 2) * norm( model_func(img_out_iter, 'notransp') - target_vector )^2;
    end
    % execution time
    execution_time_per_iter(i) = toc - start_time;
end

% done, return solution
rec_img = img_out_iter;

end

function [y] = add_nn_constraint_term_to_model_matrix(inp, flag, model_fun, lambda, beta, n, m)
% calculates the function handle for lsqr, for the augmented problem (to
% add the non-negative contraint term)
% > inp: input
% > flag: either 'transp' of 'notransp'
% > model_fun: model function
% > lambda, beta: method params
% > n, m: problem sizes (size of x and b respectively)
% returns the augmented problem function handles (see top)
% (lambda/2) * ||Au - b|| + (beta/2) * ||w - u|| = 1/2 * ||A_aug * u - b_aug||,
% where: A_aug = [sqrt(lambda) * A; sqrt(beta) * I], 
%        b_aug = [sqrt(lambda) * b; sqrt(beta) * w]

myfunc_NT = @(inp) [sqrt(lambda) * model_fun(inp, 'notransp')', sqrt(beta) * inp']'; 
myfunc_T = @(inp) [(sqrt(lambda) * model_fun(inp(1:m), 'transp'))' + sqrt(beta) * inp(m+1:end)']';

if strcmp(flag,'notransp')
    y = myfunc_NT(inp);
elseif strcmp(flag,'transp')
    y = myfunc_T(inp);
end

end

