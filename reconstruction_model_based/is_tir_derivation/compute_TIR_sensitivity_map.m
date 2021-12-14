function [SIR_sensitivity_energy_map,TIR_sensitivity_energy_map] = compute_TIR_sensitivity_map(model,plot_TIR_sensitivity_map)
%% Computes sensitivity map of the transducer based on SIR and EIR
% NOTE:
%  - plot_TIR_sensitivity_map recomputes the exact SIR for each pixel on
%    the grid and for each transducer
%%%%%%%%%%%%%%%%%%%%%%%

% Check, if SPMR is used
if ~model.includeSpmr
    warning('Sensitivity for MSOTmodel is plotted, even though SPMR is not included in definition.\n');
end

% Check for EIR
assert(model.Probe.includeEir,'No EIR included in model. Please include an EIR to allow computation of TIR sensitivity.');

% Virtual source location of phantom for all detectors
if model.useSingleSpeedOfSound
    speed_of_sound_coupling = model.Medium.speedOfSound;
else
    speed_of_sound_coupling = model.Probe.coupling.speedOfSound;
end

% Init FIELD
% Compute necessary element resolution & define elevation focused transducer
frequency_sampling_field = 5e9;            % set higher sampling for FIELD II computation
if model.useSingleSpeedOfSound
    speed_of_sound_field = model.Medium.speedOfSound;
else
    speed_of_sound_field = model.Probe.coupling.speedOfSound;
end
fov_circle_radius = sqrt(max(abs(model.Discretization.region.xPositionsFov(1,1)),abs(model.Discretization.region.xPositionsFov(1,end))).^2 + max(max(abs(model.Discretization.region.zPositionsFov(1,1)),abs(model.Discretization.region.zPositionsFov(end,1)))).^2);
distance_min_fov = ceil((model.Probe.detector.radius - fov_circle_radius)*1e3)/1e3-0.002;      % subtract offset of 2 mm from min radius
elevation_focused_transducer = init_field_II_and_get_single_toroidal_focused_transducer(model.Probe.detector.elevationRadiusOfSingleTransducer,model.Probe.detector.heightofSingleTransducer,model.Probe.detector.pitchOfTransducers,...
        model.Probe.detector.seperationBetweenTransducers,distance_min_fov,speed_of_sound_field,frequency_sampling_field);


% Zero-pad EIR
samples_zero_padded = 500;
fft_length_DAC = size(model.Probe.detector.EIR,1)+samples_zero_padded;
fft_length_field = (fft_length_DAC-1)*frequency_sampling_field/model.Probe.DAC.frequency+1;
h_EIR_hat = fft(model.Probe.detector.EIR,fft_length_DAC,1);
samples_zeros_upsampling = (fft_length_field-fft_length_DAC)/2; 
h_EIR_hat_padded = ifftshift(padarray(fftshift(h_EIR_hat),[samples_zeros_upsampling 0],0,'both'));

% Output variables (half field of view is sufficient for symetric FoV)
assert( abs(model.Discretization.region.xPositionsFov(1,1)) - abs(model.Discretization.region.xPositionsFov(1,end)) < 1e-9,'Please select a FoV with symmetry to z-axis for efficiency')
TIR_sensitivity_energy_map = zeros(size(model.Discretization.region.xPositionsFov),model.Probe.detector.numOfTransducers);
SIR_sensitivity_energy_map = zeros(size(model.Discretization.region.xPositionsFov),model.Probe.detector.numOfTransducers);

% Get active pool of parallel workers or create one (necessary for parfor-loop in fermat_principle.m)
pool = gcp;

%% Compute TIR for field of view
% Iterate over all pixel
fprintf('Current pixel index x (total number: %i, y pixels per iteration: %i): ', floor(model.Discretization.sizeOfPixelGrid(1)/2)+1,model.Discretization.sizeOfPixelGrid(2));
for index_pixel_x = 1:(floor(model.Discretization.sizeOfPixelGrid(1)/2)+1)
    fprintf('%i, ',index_pixel_x);
    
    for index_pixel_z = 1:model.Discretization.sizeOfPixelGrid(2)
   
        % Virtual source location of grid point for all detectors
        if model.useSingleSpeedOfSound
            x_phantom_virtuals = repmat(model.Discretization.region.xPositionsFov(index_pixel_z,index_pixel_x),[1,model.Probe.detector.numOfTransducers]);
            z_phantom_virtuals = repmat(model.Discretization.region.zPositionsFov(index_pixel_z,index_pixel_x),[1,model.Probe.detector.numOfTransducers]);
        else
            [x_membrane_intersections,travel_times,~] = refraction(model.Probe.membrane,model.Probe.detector.radius,model.Discretization.region.xPositionsFov(index_pixel_z,index_pixel_x),model.Discretization.region.zPositionsFov(index_pixel_z,index_pixel_x),...
                  model.Probe.detector.xPositionsOfTransducers,model.Probe.detector.zPositionsOfTransducers, model.Medium.speedOfSound,speed_of_sound_coupling);
            z_membrane_intersections = model.Probe.membrane(x_membrane_intersections); 
            x_phantom_virtuals = model.Probe.detector.xPositionsOfTransducers + (x_membrane_intersections - model.Probe.detector.xPositionsOfTransducers)./sqrt((x_membrane_intersections - model.Probe.detector.xPositionsOfTransducers).^2 + (z_membrane_intersections - model.Probe.detector.zPositionsOfTransducers).^2).*travel_times*speed_of_sound_coupling;
            z_phantom_virtuals = model.Probe.detector.zPositionsOfTransducers + (z_membrane_intersections - model.Probe.detector.zPositionsOfTransducers)./sqrt((x_membrane_intersections - model.Probe.detector.xPositionsOfTransducers).^2 + (z_membrane_intersections - model.Probe.detector.zPositionsOfTransducers).^2).*travel_times*speed_of_sound_coupling;
            nan_indices = isnan(x_membrane_intersections);

            % If xmin is not defined, set the virtual location to the actual location of the pixel
            x_phantom_virtuals(isnan(x_membrane_intersections)) = model.Discretization.region.xPositionsFov(nan_indices);
            z_phantom_virtuals(isnan(x_membrane_intersections)) = model.Discretization.region.zPositionsFov(nan_indices);
        end

        % Iterate over detectors
        for index_transducer = 1:model.Probe.detector.numOfTransducers

            % Homogeneous transformation to have source position wrt to detector
            detector_angle = ((180-model.Probe.detector.angularCoverage)/2 + (index_transducer-1)*model.Probe.detector.angularCoverage/(model.Probe.detector.numOfTransducers-1))*pi/180;
            theta_rot=detector_angle+pi/2;
            H_mat = [cos(-theta_rot) -sin(-theta_rot) 0;
                     sin(-theta_rot) cos(-theta_rot) model.Probe.detector.radius;
                     0 0 1];
            pos_virtual_rot = H_mat*[x_phantom_virtuals(index_transducer);z_phantom_virtuals(index_transducer);1];

            % SIR computation with FIELD II
            h_SIR_padded = zeros(50*frequency_sampling_field/model.Probe.DAC.frequency+1,1);
            [h_SIR,~] = calc_h(elevation_focused_transducer,[pos_virtual_rot(1) 0 pos_virtual_rot(2)]);
            h_SIR_padded(2:(length(h_SIR)+1),1) = h_SIR;
            h_SIR_padded_hat = fft(h_SIR_padded,fft_length_field,1);
            
            % Assign SIR energy (symmetry of transducers)
            SIR_sensitivity_energy_map(index_pixel_z,index_pixel_x,index_transducer) = trapz(h_SIR_padded.^2)/frequency_sampling_field;
            SIR_sensitivity_energy_map(index_pixel_z,model.Discretization.sizeOfPixelGrid(1)-index_pixel_x+1,model.Probe.detector.numOfTransducers-index_transducer+1) = SIR_sensitivity_energy_map(index_pixel_z,index_pixel_x,index_transducer);

            % Compute TIR energy (symmetry of SIR wrt to transducers)
            if size(h_EIR_hat_padded,2) > 1
                h_TIR_hat = h_SIR_padded_hat.*h_EIR_hat_padded(:,index_transducer);
                TIR_sensitivity_energy_map(index_pixel_z,index_pixel_x,index_transducer) = (trapz(abs(h_TIR_hat).^2)/fft_length_field)/frequency_sampling_field;
                h_TIR_hat = h_SIR_padded_hat.*h_EIR_hat_padded(:,model.Probe.detector.numOfTransducers-index_transducer+1);
                TIR_sensitivity_energy_map(index_pixel_z,model.Discretization.sizeOfPixelGrid(1)-index_pixel_x+1,model.Probe.detector.numOfTransducers-index_transducer+1) = (trapz(abs(h_TIR_hat).^2)/fft_length_field)/frequency_sampling_field;
            else
                h_TIR_hat = h_SIR_padded_hat.*h_EIR_hat_padded;
                TIR_sensitivity_energy_map(index_pixel_z,index_pixel_x,index_transducer) = (trapz(abs(h_TIR_hat).^2)/fft_length_field)/frequency_sampling_field;
                TIR_sensitivity_energy_map(index_pixel_z,model.Discretization.sizeOfPixelGrid(1)-index_pixel_x+1,model.Probe.detector.numOfTransducers-index_transducer+1) = TIR_sensitivity_energy_map(index_pixel_z,index_pixel_x,index_transducer);
            end

        end
    end
end

fprintf('\n');

% delete parallel pool of workers
delete(pool); 

% Plot map
if plot_TIR_sensitivity_map
    TIR_sensitivity_energy_normalized = mean(TIR_sensitivity_energy_map,3);
    TIR_sensitivity_energy_normalized = TIR_sensitivity_energy_normalized./max(TIR_sensitivity_energy_normalized(:));
    
    figure_string = ['MSOT Probe ' model.Probe.probeId ' - Normalized Sensitivity Map'];
    figure('name',figure_string);
    set(gcf, 'Units', 'centimeters', 'Position', [0 1 15 12]);
    hold on;
            
    % Sensitivity map
    imagesc(model.Discretization.region.xPositionsFov(1,:),model.Discretization.region.zPositionsFov(:,1),TIR_sensitivity_energy_normalized);
            
    % -3dB contour
    [~,c_TIR] = contour(model.Discretization.region.xPositionsFov,model.Discretization.region.zPositionsFov,TIR_sensitivity_energy_normalized,[0.7071 0.7071]);
            
    % Membrane
    if isa(model.Probe.membrane,'function_handle')
        plot(model.Discretization.region.xPositionsFov(1,:),model.Probe.membrane(model.Discretization.region.xPositionsFov(1,:)),'Color',[0, 0.4470, 0.7410]);
    end
    
    % Focus point
    p_focus_point = plot(0,0,'-xk','MarkerSize',7);
            
    hold off;
    pbapect([1 1 1]);
    xlim([model.Discretization.region.xPositionsFov(1,1) model.Discretization.region.xPositionsFov(1,end)]);
    ylim([model.Discretization.region.zPositionsFov(1,1) model.Discretization.region.zPositionsFov(1,end)]);
    view([0 -90]);
%     colormap gray;
    cb = colorbar;
    cb.TickLabelInterpreter = 'latex';
    caxis([0 1]);
    c_TIR.LineColor = 'k';%[0.8500, 0.3250, 0.0980];
    ylabel(cb, 'normalized TIR energy','Interpreter','latex', 'FontSize',14);
    set(gca,'TickLabelInterpreter', 'latex', 'FontSize',14);
    xlabel('pixel position $x_{fov}$','Interpreter','Latex', 'FontSize',14);
    ylabel('pixel position $z_{fov}$','Interpreter','Latex', 'FontSize',14);
end

end