function [rec_img, objective_per_iter, execution_time_per_iter] = nnls_sparsa_based_reg_matrix_less(model_func, Psi, Phi, target_vector, initial_img, lambda, beta, tau, nIter, nIter_inner, calculate_objective_per_iter)
% performing NNLS with some regularization function (TV or Shearlet)
% uses method from: "Total-Variation Regularization With Bound Constraints"
% Note: this method was originally derived by the authors for the TV norm,
% but we adapt it and use it also for Shearlet.
% Inputs:
% > model_func:    function for computing the forward model
% > Psi:           Psi function for SPARSA (see SPARSA documentation for details)
% > Phi:           Phi function for SPARSA (regularization function)
% > target_vector: vector of transducer data (we want the forward of rec_img to match them)
% > initial_img:   initial image - guess of the solution (usually set as vector 0)
%                  the x0 of the optimization method
% > lambda:        lambda param of the method
% > beta:          beta param of the method
% > tau:           tau regularization param of the method
% > nIter:         max. number of outer iteration
% > nIter_inner:   max. number of inner iteration (usually set below 5,
%                  depending on the problem size)
% > calculate_objective_per_iter:   Flag if objective should be
%                                   calculated at the end of every outer
%                                   iteration (one additional model evaluation)
% Outputs:
% > rec_img:                 solution - reconstructed image
% > objective_per_iter:      vector containing the objective value f(x) = ||Ax-b||^2 + tau*norm(u)
%                            after each outer iteration
% > execution_time_per_iter: vector containing the execution time after
%                            each outer iteration (in sec)
% terminates if nIter is reached or epsStepRel condition is met

% initialization
img_out_iter = zeros(size(initial_img));
lagrange_multip_for_convergence = img_out_iter;
num_target_values = length(target_vector);

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
linear_operator_augmented = @(inp) [sqrt(lambda) * model_func(inp, 'notransp')', sqrt(beta) * inp']';
% for sparsa, we need also the transpose, A_aug', which is constructed accordingly
linear_operator_augmented_T = @(inp) [(sqrt(lambda) * model_func(inp(1:num_target_values), 'transp'))' + sqrt(beta) * inp(num_target_values+1:end)']';

% outter iterations
tic;
for i = 1:nIter
    start_time = toc;
    % print interations msg
    if mod(i, 100) == 0
        disp(['iteration: ', num2str(i), ' / ', num2str(nIter)]);
    end
    % step 1: solve augmented problem with SPARSA
    % we use the Psi and Phi functions provided (see SPARSA Documentation)
    
    % apply SPARSA to get u update, which is the inner iteration image -
    % solution:
    [img_inner_iter, ~, ~, ~, ~, ~] = SpaRSA([sqrt(lambda)*target_vector', sqrt(beta)*(img_out_iter+lagrange_multip_for_convergence)']',linear_operator_augmented, tau, 'PHI', Phi, 'PSI', Psi, 'AT', linear_operator_augmented_T, 'Initialization', img_out_iter, 'MaxiterA', nIter_inner, 'StopCriterion',1, 'ToleranceA', 1e-8, 'Verbose', 0);
    
    % step 2: solve 10 from paper (w update, this is simply clipping)
    % this is the image - solution of the outer iteration
    img_out_iter = max(img_inner_iter - lagrange_multip_for_convergence, 0);
    
    % step 3: v update (see paper)
    % this is a Lagrangian multiplier ensuring convergence
    lagrange_multip_for_convergence = lagrange_multip_for_convergence + img_out_iter - img_inner_iter;
    
    % iteration end, store values
    if calculate_objective_per_iter
        objective_per_iter(i) = (lambda / 2) * norm( model_func(img_out_iter, 'notransp') - target_vector )^2 + tau * Phi(img_out_iter);
    end
    
    % execution time
    execution_time_per_iter(i) = toc - start_time;
end

% done, return solution
rec_img = img_out_iter;

end
