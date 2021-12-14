% template for linear map as required by Matlab

% given two handles A and AT, use as follows:
% afun = @(x,flag) lin_eval(x,flag,A,AT);

function val = combine_forward_and_transpose_with_reg(x, flag, A, AT, lambda, l2Reg, l2RegT, numL2RegElems)

switch flag
    case 'notransp'
        val = [vec(A(x)); lambda.*l2Reg(x)];
        val = val(:);
    case 'transp'
        val = vec(AT( x(1:(end-numL2RegElems)) )) + lambda.*l2RegT( x((end-numL2RegElems+1) : end) );
end
