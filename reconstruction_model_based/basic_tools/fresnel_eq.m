function reflection_coeff_array = fresnel_eq(incident_angle_array,speed_of_sound_below_membrane,speed_of_sound_above_membrane,density_below_membrane,density_above_membrane)

%% Calculation of the reflection coefficient via Fresnel's equation

% Fresnel's equation
impedance_above_membrane = density_above_membrane*speed_of_sound_above_membrane*cos(incident_angle_array);
impedance_below_membrane = density_below_membrane*sqrt(speed_of_sound_below_membrane^2-speed_of_sound_above_membrane^2*sin(incident_angle_array).^2);

reflection_coeff_array = ((impedance_above_membrane-impedance_below_membrane)./(impedance_above_membrane+impedance_below_membrane)).^2;