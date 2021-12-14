function [distance_stepsize_fov,distance_min_fov,distance_max_fov,angle_azimuth_stepsize_fov,angle_azimuth_max_fov,elevation_stepsize_fov,elevation_max_fov,single_element_resolution_FieldII_min] = ...
    calculate_SIR_params(detector_array_radius,x_fov_min,x_fov_max,...
                                      y_fov_min,y_fov_max,z_fov_max,desired_resolution_fov,speed_of_sound_coupling,frequency_sampling_field)

%% Calculation of polar field of view, single element resolution and SIR length 
%% of a single MSOT detector element for an environment with single speed of sound
% NOTE: 
%  - polarSIR dimensions approximated by circle through edges of field of view in plane
%  - Snap to the grid in models causes ripple artefacts in sinogram, if SIR computed on grid which is too sparse
%  - Calculation of element resolution is based on following paper:
%    J.A. Jensen, N.B. Svendsen. "Calculation of Pressure Fields from Arbitrarily Shaped, Apodized, and Excited Ultrasound Transducers." 
%    IEEE Transactions on Ultrasonics, Ferroelectrics and Frequency Control, vol. 39, no. 2, pp. 262-267, 1991
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate radius and azimuthal angle boundaries to match parameters 
% Principle: draw circle through FoV edge (largest distance from focus), centered in focus
fov_circle_radius = sqrt(max(abs(x_fov_min),abs(x_fov_max)).^2 + max(abs(y_fov_min),abs(y_fov_max)).^2);
distance_min_fov = ceil((detector_array_radius - fov_circle_radius)*1e3)/1e3-0.002;      % subtract offset of 2 mm from min radius
distance_max_fov = ceil((detector_array_radius + fov_circle_radius)*1e3)/1e3+0.002;      % add offset of 2 mm to max radius
angle_azimuth_max_fov = ceil(asin(fov_circle_radius/detector_array_radius)/pi*180)+2;     % add offset of 2 deg to max angle 
elevation_max_fov = z_fov_max;

% Calculate step sizes to match resolution
distance_stepsize_fov = desired_resolution_fov/4;                                                 % min. 4x desired resolution of grid in FoV required in r-dimension 
angle_azimuth_stepsize_fov = ceil(atan2(desired_resolution_fov/4,distance_max_fov)/pi*180*100)/100;
elevation_stepsize_fov = desired_resolution_fov/2;                                                 % 2x desired resolution of grid in FoV required in z-dimension

% Calculate length of discretized square in detector 
% (according to eq. (17) in paper, all points of FoV have to be in far-field region of detector element)
single_element_resolution_FieldII_min = 50e-6;
elem_res_min = sqrt(4*distance_min_fov*speed_of_sound_coupling/frequency_sampling_field)/10;
if (elem_res_min < single_element_resolution_FieldII_min)
    single_element_resolution_FieldII_min = round(elem_res_min*1e6)/1e6;
end

end