function travel_time_array = calculate_travel_time_with_dual_speed_of_sound_without_intersec(membrane_model,x_fov_array,z_fov_array,x_detector_array,z_detector_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane)

%% Calculation of time of flight for a dual SoS model, using the C-codes of the Fast Marching Implementation by Dirk-Jan Kroon (https://de.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching)

%% preliminaries
% conversion from detector coordinate to indices [ix,iy]
% (this introduces a small error, depending on the discretization)
delta_x = abs(x_fov_array(1,2)-x_fov_array(1,1));
delta_z = abs(z_fov_array(2,1)-z_fov_array(1,1));

x_offset = min(x_fov_array(:));
z_offset = min(z_fov_array(:));

index_x = round((x_detector_array-x_offset)/delta_x)+1;
index_z = round((z_detector_array-z_offset)/delta_z)+1;

% setup of the speed of sound distribution
speed_of_sound_distrib = speed_of_sound_above_membrane*(z_fov_array<=membrane_model(x_fov_array)) + speed_of_sound_below_membrane*(z_fov_array>membrane_model(x_fov_array));

% initialization
travel_time_array = zeros([size(x_fov_array),length(x_detector_array)]);

%% time of flight calculation, loop over detectors
for k = 1:length(x_detector_array)
    % fast marching, using second derivatives, and cross neighbors
    travel_time_array(:,:,k) = delta_x*msfm2d(speed_of_sound_distrib, [index_z(k),index_x(k)]',true,true);
end