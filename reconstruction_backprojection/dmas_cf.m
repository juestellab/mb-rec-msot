% DMAS_CF: Delay-Multiply-And-Sum with Coherence Factor
%
% [1] Matrone G, Savoia AS, Caliano G, Magenes G.
% The delay multiply and sum beamforming algorithm in ultrasound B-mode medical imaging.
% % IEEE transactions on medical imaging, 34(4):940-949, 2014.
%
% [2] Park J, Jeon S, Meng J, Song L, Lee JS, Kim C.
% Delay-multiply-and-sum-based synthetic aperture focusing in photoacoustic microscopy.
% Journal of biomedical optics, 21(3):036010, 2016.
%
% [3] Mozaffarzadeh M, Mahloojifar A, Orooji M, Adabi S, Nasiriavanaki M.
% Double-stage delay multiply and sum beamforming algorithm: Application to linear-array photoacoustic imaging.
% IEEE Transactions on Biomedical Engineering, 65(1):31-42, 2017.
%
% [4] Mozaffarzadeh M, Sadeghi M, Mahloojifar A, Orooji M.
% Double-stage delay multiply and sum beamforming algorithm applied to ultrasound medical imaging.
% Ultrasound in Medicine & Biology, 44(3):677-86, 2018.

function [p0_rec_dmas_cf, p0_rec_dmas, cf, p0_rec_das] = dmas_cf(...
  sigMat,...
  fov_x,...
  fov_z,...
  speed_of_sound,...
  sampling_frequency,...
  cropped_or_unrecorded_signals_at_sinogram_start,...
  angular_coverage,...
  detector_radius)

  % Derive required properties from input parameters
  dt = 1/sampling_frequency;
  number_of_samples = size(sigMat, 1);
  number_of_transducers = size(sigMat, 2);
  angular_offset = (180 - angular_coverage) / 2;
  transducer_angles = -(angular_offset : angular_coverage/(number_of_transducers-1) : (180-angular_offset) ) * pi/180;
  aquisition_times_in_samples = (1:number_of_samples) + cropped_or_unrecorded_signals_at_sinogram_start;
  aquisition_times_in_seconds = aquisition_times_in_samples * dt;

  sigMat = sigMat / max(abs(sigMat(:)));
  
  % Take the signed square root of the absolute sigMat values, to preserve dimensionality after pairwise multiplication
  % "\hat{k} = sign(k) * sqrt(abs(k))"
  sigMat_sqrt = sign(sigMat) .* sqrt(abs(sigMat));

  % Initialisation of the final outputs
  p0_rec_dmas_cf = zeros(size(fov_x,1), size(fov_x,2), size(sigMat,3));
  p0_rec_dmas = zeros(size(fov_x,1), size(fov_x,2), size(sigMat,3));
  p0_rec_das = zeros(size(fov_x,1), size(fov_x,2), size(sigMat,3));
  cf = zeros(size(fov_x,1), size(fov_x,2), size(sigMat,3));

  
  for wavelength = 1:size(sigMat,3)

    % DMAS involves a series of pairwise multiplications:
    % "dmas(k) = \sum_{i=1}^{i=M-1} \sum{j=i+1}^{j=M} \hat{k}_i * \hat{k}_j"
    % However this operation can be simplified as follows:
    % "dmas(k) = 0.5 * (sum(k)^2 - sum(k^2))" --> We will store "sum_of_squares" and "squared_sum"
    p0_rec_sum_of_squares = zeros(size(fov_x));
    p0_rec_squared_sum = zeros(size(fov_x));

    sigMat_wl = sigMat(:,:,wavelength);
    sigMat_wl_sqrt = sigMat_sqrt(:,:,wavelength);    
  
    % Coherence factor (Note that the numerator is calculated exactly like a DAS operation)
    cf_numerator = zeros(size(fov_x,1), size(fov_x,2));
    cf_denominator = zeros(size(fov_x,1), size(fov_x,2));    
    
    for p = 1:number_of_transducers
      % Calculate time needed for propagation from every pixel to detector and convert to indices in the time array (nearest interpolation)
      T = round((...
        sqrt(abs(fov_x(:)-cos(transducer_angles(p))*detector_radius).^2 + abs(fov_z(:)-sin(transducer_angles(p))*detector_radius).^2)...
        / speed_of_sound...
        - aquisition_times_in_seconds(1)...
        )*sampling_frequency + 1 ...
        );
  
      % Add all values in sigMat as indexed by T
      % The terms corresponding to "sum_of_squares" must he squared here
      % The terms corresponding to "squared_sum" will be squared once the sum is over (after the FOR loop)      
      p0_rec_sum_of_squares(:) = p0_rec_sum_of_squares(:) + sigMat_wl_sqrt(T,p) .^2; % Squaring
      p0_rec_squared_sum(:) = p0_rec_squared_sum(:) + sigMat_wl_sqrt(T,p); % Will be squared after the FOR loop
      
      cf_numerator(:) = cf_numerator(:) + sigMat_wl(T,p); % Will be squared after the FOR loop (NOTE: This is also the process for the simple DAS)
      cf_denominator(:) = cf_denominator(:) + sigMat_wl(T,p) .^2; % Squaring
      
    end

    % Simple DAS
    p0_rec_das(:,:,wavelength) = cf_numerator;
    
    % Square the sums
    p0_rec_squared_sum  = p0_rec_squared_sum .^2;
    cf_numerator = cf_numerator .^2;
    
    % We build the DMAS result s.t.: "dmas(y) = 0.5 * (sum(y)^2 - sum(y^2))"
    p0_rec_dmas(:,:,wavelength) = 0.5 * (p0_rec_squared_sum - p0_rec_sum_of_squares);

    % We build the CF, which involves computing the pixel-wise ratio of a squared sum divided by a sum of squares
    % The underlying principle of CF is analogous to the Cauchy-Schwarz inequalit
    % "CF(k) = sum(k)^2 / (M * sum(k^2))"
    cf(:,:,wavelength) = cf_numerator ./ (1e-6 + number_of_transducers * cf_denominator);
   
  end
  
  % We apply the CF to the DMAS result
  p0_rec_dmas_cf = p0_rec_dmas .* cf;
  
end
