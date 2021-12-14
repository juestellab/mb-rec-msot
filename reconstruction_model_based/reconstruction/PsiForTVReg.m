function [y] = PsiForTVReg(x, tau,B,Bt,BtB)
% apply the Psi function on all wavelengths
y = zeros(size(x,1)*size(x,2), size(x,3));

parfor i = 1:size(x,3)
    y(:,i) = SB_ITV2(x(:,:,i),tau,B,Bt,BtB);
end
y = y(:);
end