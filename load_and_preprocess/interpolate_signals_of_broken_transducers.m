function sigMat_interpolated = interpolate_signals_of_broken_transducers(sigMat, broken_transducers)

if isempty(broken_transducers)
    sigMat_interpolated=sigMat;
else
    non_broken_transducers = setdiff(1:size(sigMat,2), broken_transducers);

    [T, det] = meshgrid(1:size(sigMat,1), 1:size(sigMat,2));
    [T_non_broken, det_non_broken] = meshgrid(1:size(sigMat,1), non_broken_transducers);


    sigMat_interpolated = zeros(size(sigMat));
    for f=1:size(sigMat,4)
        for wl=1:size(sigMat,3)
            sigMat_interpolated(:,:,wl,f) = interp2(T_non_broken, det_non_broken, sigMat(:,non_broken_transducers,wl, f)', T, det, 'spline')';
        end
    end

end
end