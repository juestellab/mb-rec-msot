function [travel_time_array,transmission_coeff_array,pos_fov_virtual_array_struct] = calculate_SoS_model(membrane_model,x_fov_array,y_fov_array,x_detector_array,y_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane,density_below_membrane,density_above_membrane)

%% Calculation of the travel times and transmission coefficients for refraction at an interface
% Restrictions: values of the transmission coefficient are inexact very
%               close to the membrane due to numerical errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resize mesh-grid
x_fov_array_size = size(x_fov_array);
x_fov_array = x_fov_array(:);
y_fov_array = y_fov_array(:);

travel_time_array = zeros([size(x_fov_array),length(x_detector_array(:))]);
transmission_coeff_array = zeros(size(travel_time_array));
x_fov_virtual_array = zeros(size(travel_time_array));
y_fov_virtual_array = zeros(size(travel_time_array));

% With for-loop: execution time 234.512617 s (on grid with dim: 41 x 41 x 256)
% With parfor-loop: execution time 72.622022 s (on grid with dim: 41 x 41 x 256)
parfor j = 1:length(x_fov_array)
% for j = 1:length(x_resh)      % Only for debugging purposes
        [x_membrane_intersection_array,travel_time_array(j,:),~,transmission_coeff_array(j,:)] = refraction(membrane_model,x_fov_array(j),y_fov_array(j),x_detector_array,y_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane,density_below_membrane,density_above_membrane);
        if sum(isnan(x_membrane_intersection_array)) >= 1
            x_fov_virtual_array(j,:) = ones(length(x_detector_array),1)*x_fov_array(j);
            y_fov_virtual_array(j,:) = ones(length(x_detector_array),1)*y_fov_array(j);
        else
            delta_x_array = x_detector_array-x_membrane_intersection_array;
            delta_y_array = y_detector_array-membrane_model(x_membrane_intersection_array);
            delta_distance_array = sqrt(delta_x_array.^2+delta_y_array.^2);
            x_fov_virtual_array(j,:) = x_detector_array - delta_x_array./delta_distance_array.*travel_time_array(j,:)*speed_of_sound_above_membrane;
            y_fov_virtual_array(j,:) = y_detector_array - delta_y_array./delta_distance_array.*travel_time_array(j,:)*speed_of_sound_above_membrane;
        end
end

% reshape matrices
travel_time_array = reshape(travel_time_array,[x_fov_array_size,length(x_detector_array)]);
transmission_coeff_array = reshape(transmission_coeff_array,[x_fov_array_size,length(x_detector_array)]);
pos_fov_virtual_array_struct.x = reshape(x_fov_virtual_array,[x_fov_array_size,length(x_detector_array)]);
pos_fov_virtual_array_struct.y = reshape(y_fov_virtual_array,[x_fov_array_size,length(x_detector_array)]);