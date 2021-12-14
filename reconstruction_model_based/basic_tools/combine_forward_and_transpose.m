% template for linear map as required by Matlab

% given two handles A and AT, use as follows:
% afun = @(x,flag) lin_eval(x,flag,A,AT);

function val = combine_forward_and_transpose(x,flag,A,AT)

switch flag
    case 'notransp'
        val = A(x);
        val = val(:);
    case 'transp'
        val = AT(x);
        val = val(:);
end
