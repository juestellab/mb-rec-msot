function K = bp_kernel(t)

Nt = length(t);
dlambda = pi/(t(end)-t(1));
lambda = reshape(dlambda:dlambda:length(t)*dlambda,[1,1,Nt]);
H0 = besselh(0,bsxfun(@times,t',lambda));

K = zeros(Nt);
for j = 1:Nt
    K(j,:) = trapz(squeeze(lambda),H0(j,:,:).*conj(H0).*lambda,3);
end
K = -imag(K);