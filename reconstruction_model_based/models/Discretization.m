classdef Discretization
    properties
        region
        pixelModel
        sizeOfPixelGrid
        numberOfPixels
        prefilter
        prefilterConj
    end
    
    methods (Access=public)
        function discretization = Discretization(field_of_view, size_of_pixel_grid, custom_pixel_model)
            if nargin < 3
                custom_pixel_model = [];
            elseif ~isa(custom_pixel_model,'function_handle')
                error('Custom pixel model must be function handle.');
            end
            
            % Pixel grid
            if length(size_of_pixel_grid) > 1
                discretization.sizeOfPixelGrid = size_of_pixel_grid;
            else
                discretization.sizeOfPixelGrid = [size_of_pixel_grid, size_of_pixel_grid];  % Discr.sizeOfPixelGrid(1): N_x, Discr.sizeOfPixelGrid(2): N_z
            end
            
            % Field of view
            if (length(field_of_view) <= 2)
                x_vec_fov = linspace(field_of_view(1), field_of_view(2), discretization.sizeOfPixelGrid(1));
                z_vec_fov = linspace(field_of_view(1), field_of_view(2), discretization.sizeOfPixelGrid(2));
            else
                x_vec_fov =  linspace(field_of_view(1), field_of_view(2), discretization.sizeOfPixelGrid(1));
                z_vec_fov =  linspace(field_of_view(3), field_of_view(4), discretization.sizeOfPixelGrid(2));
            end
            [discretization.region.xPositionsFov,discretization.region.zPositionsFov] = meshgrid(x_vec_fov,z_vec_fov);
            discretization.numberOfPixels = length(discretization.region.xPositionsFov(:));
            
            % Pixel model
            if ~isempty(custom_pixel_model)
                assert(isa(custom_pixel_model,'function_handle'),'custom_pixel_model has to be a function handle');
                discretization.pixelModel = custom_pixel_model;
            else
                resolution_fov = min((discretization.region.xPositionsFov(1,2)-discretization.region.xPositionsFov(1,1)),(discretization.region.zPositionsFov(2,1)-discretization.region.zPositionsFov(1,1)));
                discretization.pixelModel = @(r) gaussian_pixel_model(r,round(resolution_fov/2*1e6)/1e6);
            end
            
            % Prefilter for projection into shift-invariant space
            pixel_on_grid = discretization.evaluate_pixel_model_on_grid(1);
            [~,~,prefilter] = shift_inv_proj(zeros(discretization.sizeOfPixelGrid(2), discretization.sizeOfPixelGrid(1)), [1 1], fftshift(pixel_on_grid));
            
            discretization.prefilter = prefilter;
            discretization.prefilterConj =  conj(prefilter);
        end
        
        function coeffs = project_into_shift_invariant_space(discretization, initial_pressure)
            coeffs = reshape(initial_pressure, discretization.sizeOfPixelGrid(2), discretization.sizeOfPixelGrid(1),[]);
            
            % This filtering performs a circular instead of a linear convolution (no zero-padding of signal).
            % However, here we consider the introduced errors towards the image boundaries as negligible.
            for i=1:size(coeffs,3)
                coeffs(:,:,i) = ifftn(fftn(coeffs(:,:,i)).*discretization.prefilter,'symmetric');
            end
        end
        
        function coeffs_transp = transpose_project_into_shift_invariant_space(discretization, initial_pressure_transp)
            coeffs_transp = reshape(initial_pressure_transp, discretization.sizeOfPixelGrid(2), discretization.sizeOfPixelGrid(1),[]);
            
            % This filtering performs a circular instead of a linear convolution (no zero-padding of signal).
            % However, here we consider the introduced errors towards the image boundaries as negligible.
            for i=1:size(coeffs_transp,3)
                coeffs_transp(:,:,i) = ifftn(fftn(coeffs_transp(:,:,i)).*discretization.prefilterConj,'symmetric');
            end
        end
        
        function pixel_on_grid = evaluate_pixel_model_on_grid(discretization, refinement_factor)
            grid_N = refinement_factor .* (discretization.sizeOfPixelGrid-1) + 1;
            grid_spacing = abs([discretization.region.xPositionsFov(1,2)-discretization.region.xPositionsFov(1,1), discretization.region.zPositionsFov(2,1)-discretization.region.zPositionsFov(1,1)]) / refinement_factor;
            
            % The grid, on which the pixel model is evaluated, must be centered
            % so that the peaks of fftshift(pixel) is located at [1,1]
            linspace_for_pixel_x = grid_spacing(1) .* (floor(-grid_N(1)/2) : 1 : floor(grid_N(1)/2-1));
            linspace_for_pixel_z = grid_spacing(2) .* (floor(-grid_N(2)/2) : 1 : floor(grid_N(2)/2-1));
            
            [x_grid_for_pixel, z_grid_for_pixel] = meshgrid(linspace_for_pixel_x, linspace_for_pixel_z);
            r = sqrt(x_grid_for_pixel.^2+z_grid_for_pixel.^2);
            pixel_on_grid = discretization.pixelModel(r);
        end
        
        function image_refined = evaluate_image_on_refined_grid(discretization, image, refinement_factor, clip_negative_values)
            % Get image coeffs of on original grid
            coeffs = discretization.project_into_shift_invariant_space(image);
            
            % Evaluate pixel model on refined grid
            pixel_on_refined_grid = discretization.evaluate_pixel_model_on_grid(refinement_factor);
            
            % Copy image coeffs to refined grid
            coeffs_refined_grid = zeros(size(pixel_on_refined_grid, 1), size(pixel_on_refined_grid, 2), size(image, 3));
            coeffs_refined_grid(1:refinement_factor:end, 1:refinement_factor:end, :) = coeffs;
            
            % Evaulate image on refined grid
            % This filtering performs a circular instead of a linear convolution (no zero-padding of signal).
            % However, here we consider the introduced errors towards the image boundaries as negligible.
            image_refined = ifft2(fft2(coeffs_refined_grid) .* fft2(fftshift(pixel_on_refined_grid)),'symmetric');
            
            if clip_negative_values
                image_refined = max(0, image_refined);
            end
        end
        
    end    
end

% Gaussian pixel model
function y = gaussian_pixel_model(r,standard_deviation)
    if nargin < 2
        error('Standard deviation for Gaussian pixel model not set.');
    end
    y = exp(-r.^2/(2*standard_deviation.^2));
end
