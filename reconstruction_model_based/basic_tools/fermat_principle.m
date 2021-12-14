function [x_membrane_intersection_array,travel_time_array] = fermat_principle(membrane_model,radius_detector_array,x_fov,z_fov,x_detector_array,z_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane)

%% Calculation of shortest path and intersection point for acoustic refraction at an interface
% Restrictions: only gives reasonable results, if [x_fov,y_fov] and [x_detector_array,y_detector_array] are on different sides of m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
x_membrane_intersection_array = zeros(length(x_fov(:)),length(x_detector_array(:)));
travel_time_array = zeros(length(x_fov(:)),length(x_detector_array(:)));

% minimization of travel time for every pair of points [x_fov,z_fov], [x_detect,z_detect]
for j = 1:length(x_fov(:))
    % Get single values (avoid broadcast variables)
    x_fov_j = x_fov(j);
    z_fov_j = z_fov(j);
    
    parfor k = 1:length(x_detector_array(:))
        
        % time needed to travel from [x_fov,z_fov] to [x_detect,z_detect] with given speeds
        % of sound sound_speed_above_membrane and sound_speed_below_membrane
        travel_time = @(x_mem) norm([x_fov_j,z_fov_j]-[x_mem,membrane_model(x_mem)],2)/speed_of_sound_below_membrane + norm([x_mem,membrane_model(x_mem)]-[x_detector_array(k),z_detector_array(k)],2)/speed_of_sound_above_membrane;
        
        % minimization of travel time over the interval [-radius_detector_array,radius_detector_array]
        [x_mem_j(k),travel_time_min_j(k)] = fminbnd(travel_time,-radius_detector_array,radius_detector_array);
        
        
    end
    
    % Save variables
    x_membrane_intersection_array(j,:) = x_mem_j;
    travel_time_array(j,:) = travel_time_min_j;
end