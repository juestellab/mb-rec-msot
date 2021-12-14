function [x_membrane_intersection_array,travel_time_array,incident_angle_array] = refraction(membrane_model,radius_detector_array,x_fov,z_fov,x_detector_array,z_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane)

%% Simulation of acoustic refraction at a simple interface
% Restrictions: values of the transmission coefficient are inexact very
%               close to the membrane due to numerical errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if z_fov > membrane_model(x_fov) % if the point is on the other side, i.e. below the interface
    
    % calculation of shortest path via Fermat's principle (minimize travel time)
    [x_membrane_intersection_array,travel_time_array] = fermat_principle(membrane_model,radius_detector_array,x_fov,z_fov,x_detector_array,z_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane);
    
    % drop a perpendicular at the point of intersection
    z_membrane_intersection_array = membrane_model(x_membrane_intersection_array);
    [x_vec_perp_membrane,z_vec_perp_membrane] = drop_vector_perpend_membrane(membrane_model,x_membrane_intersection_array,0.0001);
    
    % calculate incident angle
    incident_angle_array = acos((x_vec_perp_membrane.*(x_fov-x_membrane_intersection_array) + z_vec_perp_membrane.*(z_fov-z_membrane_intersection_array))./sqrt(x_vec_perp_membrane.^2+z_vec_perp_membrane.^2)./sqrt((x_fov-x_membrane_intersection_array).^2+(z_fov-z_membrane_intersection_array).^2));

else % if the point is on the same side, i.e. above the interface
    
    x_membrane_intersection_array = NaN(size(x_detector_array)); % no intersection point
    travel_time_array = sqrt((x_detector_array-x_fov).^2 + (z_detector_array-z_fov).^2)./speed_of_sound_above_membrane; % time for direct travel
    incident_angle_array = NaN(size(x_detector_array)); % incident angle    
end


