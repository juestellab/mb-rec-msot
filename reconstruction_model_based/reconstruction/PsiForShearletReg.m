function [y] = PsiForShearletReg(x, tau, shearletTrafo, shearletRecon)
% apply the Psi function on all wavelengths
Psi_ow = @(inp,tt) softThresholdLinearTransformation(inp,tt,shearletTrafo,shearletRecon);
y = zeros(size(x,1)*size(x,2), size(x,3));

parfor i = 1:size(x,3)
    y(:,i) = Psi_ow(x(:,:,i), tau);
end
y = y(:);
end