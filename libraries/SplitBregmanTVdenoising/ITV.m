function phi_ITV = ITV(x)

% isotropic TV norm

[dx,dy] = gradient(x);
phi_ITV = sum(sqrt(abs(dx(:)).^2+abs(dy(:)).^2));