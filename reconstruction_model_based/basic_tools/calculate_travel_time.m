function [x_membrane_intersection_array,travel_time_array] = calculate_travel_time(membrane_model,radius_detector_array,x_fov_array,z_fov_array,x_detector_array,z_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane)

%% Calculation of interface intersection points, travel times and transmission coefficients for refraction at an interface
% Restrictions: values of the transmission coefficient are inexact very
%               close to the membrane due to numerical errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resize mesh-grid
x_fov_array_size = size(x_fov_array);
x_fov_array = x_fov_array(:);
z_fov_array = z_fov_array(:);

x_membrane_intersection_array = zeros([length(x_fov_array),length(x_detector_array(:))]);
travel_time_array = zeros([length(x_fov_array),length(x_detector_array(:))]);

parfor j = 1:length(x_fov_array)
    [x_membrane_intersection_array(j,:),travel_time_array(j,:),~] = refraction(membrane_model,radius_detector_array,x_fov_array(j),z_fov_array(j),x_detector_array,z_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane);
end

travel_time_array = reshape(travel_time_array,[x_fov_array_size,length(x_detector_array)]);
x_membrane_intersection_array = reshape(x_membrane_intersection_array,[x_fov_array_size,length(x_detector_array)]);
