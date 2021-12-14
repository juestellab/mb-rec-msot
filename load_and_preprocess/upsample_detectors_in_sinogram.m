function [sinograms_upsampled] = upsample_detectors_in_sinogram(sinograms, factor)
num_upsampled_detectors = size(sinograms,2)*factor;

[T,det] = meshgrid(1:size(sinograms,1), 1:factor:num_upsampled_detectors);
[T2,det2] = meshgrid(1:size(sinograms,1), 1:num_upsampled_detectors);

sinograms_upsampled = zeros(size(sinograms,1), num_upsampled_detectors, size(sinograms,3), size(sinograms,4));
for f=1:size(sinograms,4)
    for wl=1:size(sinograms,3)
        sinograms_upsampled(:,:,wl,f) = interp2(T,det,sinograms(:,:,wl, f)',T2,det2, 'spline')';
    end
end
end

