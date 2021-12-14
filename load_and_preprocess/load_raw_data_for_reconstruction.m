function [data] = load_raw_data_for_reconstruction(model)
%LOAD_RAW_DATA_FOR_RECONSTRUCTION Function for data loading
% This function currently contains a placeholder implementation that
% simulates recorded signals with a given forward model.
%
% TODO: Load data based on input parameters

test_img = 0.01*ones(model.Discretization.sizeOfPixelGrid);
test_img(model.Discretization.region.zPositionsFov<0.015) = 0.1;
r = sqrt(model.Discretization.region.xPositionsFov.^2+model.Discretization.region.zPositionsFov.^2);
test_img(r<0.001) = 1;
test_img(10:20,23:30) = 1;
test_sig = model.Funcs.applyForward(test_img);
test_sig = cat(1, zeros(model.DataPreprocessing.numCroppedSamplesAtSinogramStart, model.Probe.detector.numOfTransducers), test_sig);
data = test_sig + 0.01 * randn(size(test_sig));
end

