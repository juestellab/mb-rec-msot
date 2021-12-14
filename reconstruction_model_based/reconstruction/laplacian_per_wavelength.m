function [x] = laplacian_per_wavelength(x)
%LAPLACIAN_PER_WAVELENGTH 

for wl =1:size(x,3)
    x(:,:,wl) = del2(x(:,:,wl),1);
end

x = x(:);
end

