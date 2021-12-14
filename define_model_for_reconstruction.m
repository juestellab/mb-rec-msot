function [model] = define_model_for_reconstruction(field_of_view, number_of_grid_points_fov, device_probe_id, use_eir, use_indiv_eir, use_sir, use_single_speed_of_sound, speed_of_sound_tissue, num_cropped_samples_at_sinogram_start, filt_cutoff_min, filt_cutoff_max, model_normalization_factor)
if nargin<12
    model_normalization_factor = [];
end

discretization = Discretization(field_of_view, number_of_grid_points_fov);
probe = Probe(device_probe_id, use_eir, use_indiv_eir);
model = MSOTmodel(probe, discretization, speed_of_sound_tissue ,use_sir, use_single_speed_of_sound, num_cropped_samples_at_sinogram_start, filt_cutoff_min, filt_cutoff_max, model_normalization_factor);

end

