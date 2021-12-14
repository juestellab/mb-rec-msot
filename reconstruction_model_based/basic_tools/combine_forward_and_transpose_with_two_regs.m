% template for linear map as required by Matlab

% given two handles A and AT, use as follows:
% afun = @(x,flag) lin_eval(x,flag,A,AT);

function val = combine_forward_and_transpose_with_two_regs(x, flag, A, AT, lambda, l2Reg, l2RegT, numL2RegElems, lambda_2, l2Reg_2, l2RegT_2)

switch flag
    case 'notransp'
        val = [vec(A(x)); lambda.*l2Reg(x); lambda_2*l2Reg_2(x)];
        val = val(:);
    case 'transp'
        val = vec(AT( x(1:(end-numL2RegElems)) )) + lambda.*l2RegT( x((end-numL2RegElems+1) : (end-(numL2RegElems/2))) ) + lambda_2.*l2RegT_2( x((end-(numL2RegElems/2)+1) : end) );
end
