function [x_vec_perp_membrane,z_vec_perp_membrane] = drop_vector_perpend_membrane(membrane_model,x_membrane_intersection,delta_x_membrane_intersection)
%% Drop a perpendicular on a line at a point
% Restrictions: 
%  - delta_x_membrane_intersection should be small value (tangent estimation)
%  - line should be smooth for good performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the perpendicular
% based on the estimated tangent via the symmetric difference quotient
x_vec_perp_membrane = membrane_model(x_membrane_intersection+delta_x_membrane_intersection) - membrane_model(x_membrane_intersection-delta_x_membrane_intersection);
z_vec_perp_membrane = 2*delta_x_membrane_intersection;

% normalization
vec_perp_membrane_norm = sqrt(x_vec_perp_membrane.^2 + z_vec_perp_membrane.^2);

x_vec_perp_membrane = x_vec_perp_membrane./vec_perp_membrane_norm;
z_vec_perp_membrane = z_vec_perp_membrane./vec_perp_membrane_norm;