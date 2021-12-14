function sinogram_mask = get_sinogram_mask_of_model_reach(model)

img = ones(model.Discretization.sizeOfPixelGrid);

forward = model.Funcs.applyForward(img);
sinogram_mask = (abs(forward)> 1e-10);

end