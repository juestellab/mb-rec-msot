function elevation_focused_transducer = init_field_II_and_get_single_toroidal_focused_transducer(radius_elevation_transducer,height_transducer,pitch_of_transducers,seperation_between_transducers,distance_fov_min,speed_of_sound_coupling,frequency_sampling_field)
%% Initialize FIELD II for SIR computation, define aperture of single transducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if FIELD II is somewhere in MATLAB path
assert(~isempty(which('field_init.m')),...
	'FIELD II is either not in the MATLAB path or installed. SIR cannot be computed without FIELD II.');

% Init FIELD II, set params & get transducer model
field_init(0);
evalc('set_field(''fs'',frequency_sampling_field)');    % use evalc functionality to catch unnecessary warning printout by FIELD II that pulses have to be reset
set_field('use_triangles',0);
set_field('use_att',0);
set_field('c',speed_of_sound_coupling);

% Set transducer element parameters
number_of_transducers = 1; % Number of elements set to 1 as only a single focused transducer is computed
electronic_focus = [0 0 0]/1000; % initial electronic focus (no focusing by phasing of signals applied)
    
% Compute maximum size of discretized elements for far field approximation to hold
size_transducer_discretization_max = sqrt(4*distance_fov_min*speed_of_sound_coupling/frequency_sampling_field)/10;
if size_transducer_discretization_max > 5e-6
    size_transducer_discretization = size_transducer_discretization_max - mod(size_transducer_discretization_max,5e-6);
else
    %% TODO: UPDATE THIS TO HOLD FOR ALL DISTANCE_FOV_MIN
    size_transducer_discretization = round(size_transducer_discretization_max,2,'significant');
end
    
% Discretize transducer
element_number_azimuth = ceil((pitch_of_transducers-seperation_between_transducers)/size_transducer_discretization); % x corresponds to direction of width in FIELD II
if mod(element_number_azimuth,2)==0
  	element_number_azimuth=element_number_azimuth+1;
end
element_number_elevation = ceil(height_transducer/size_transducer_discretization); % y corresponds to direction of height in FIELD II
if mod(element_number_elevation,2)==0
 	element_number_elevation=element_number_elevation+1;
end

% Pointer on single detector aperture
elevation_focused_transducer = xdc_focused_array(number_of_transducers,(pitch_of_transducers-seperation_between_transducers), ...
                            height_transducer,seperation_between_transducers,radius_elevation_transducer,...
                            element_number_azimuth, element_number_elevation, electronic_focus);
end
