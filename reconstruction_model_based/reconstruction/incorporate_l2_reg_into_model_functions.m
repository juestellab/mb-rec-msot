function [applyModel, targetVector] = incorporate_l2_reg_into_model_functions(data_filt, N, forward, transpose, lambdaL2, RegL2, RegL2T, lambdaL2_2, RegL2_2, RegL2T_2)
%incorporate_l2_reg_into_model_functions Incorporate l2 regularization in
%the an apply model function and a corresponding target vector.

addL2Reg = exist('lambdaL2', 'var') && ~isempty(lambdaL2)  && lambdaL2 > 0.0;
addL2Reg_2 = exist('lambdaL2_2', 'var') && ~isempty(lambdaL2_2) && (lambdaL2_2 > 0.0);

if addL2Reg && addL2Reg_2
    numElemL2Reg = 2*N(2)*N(1)*size(data_filt,3);
    applyModel = @(x,flag) combine_forward_and_transpose_with_two_regs(x, flag, forward, transpose, sqrt(lambdaL2), RegL2, RegL2T, numElemL2Reg, sqrt(lambdaL2_2), RegL2_2, RegL2T_2);

    targetVector = [data_filt(:); zeros(numElemL2Reg,1)];

elseif addL2Reg
    numElemL2Reg = N(2)*N(1)*size(data_filt,3);
    applyModel = @(x,flag) combine_forward_and_transpose_with_reg(x, flag, forward, transpose, sqrt(lambdaL2), RegL2, RegL2T, numElemL2Reg);
    targetVector = [data_filt(:); zeros(numElemL2Reg,1)];

else
    applyModel = @(x,flag) combine_forward_and_transpose(x, flag, forward, transpose);
    targetVector = data_filt(:);
end
end

