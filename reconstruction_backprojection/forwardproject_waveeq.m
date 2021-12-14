function sigMat = forwardproject_waveeq(p0,X,Y,c0,dev,psf)

% preliminaries, initialization
a = 1;
dt = 1/dev.fsample/a;
Nt = a*(dev.samples + dev.delayNumber);
t = dt:dt:Nt*dt;
dR = sqrt((X(1,2) - X(1,1))^2 + (Y(2)-Y(1))^2);

sigMat = zeros(Nt,dev.proj);

% loop over transducers
for p = 1:dev.proj
    
    % Calculate time needed for propagation from every pixel to detector and convert to indices in the time array (nearest interpolation)
    T = round((sqrt(abs(X(:)-cos(dev.angle_sensor(p))*dev.r_sensor).^2 + abs(Y(:)-sin(dev.angle_sensor(p))*dev.r_sensor).^2)/c0 -t(1))*a*dev.fsample + 1);
    % handle exceptions
    T(T<1) = 1; T(T>Nt) = Nt;
    
    % mean of all values in Recon with equal index in T
    sigMat(1:max(T),p) = accumarray(T,p0(:))*dR;
    
end

% postfiltering
sigMat = convn(sigMat,psf,'same')/sum(psf(:));

% conversion of radial integrals to pressure
sigMat = c0*gradient(sigMat'.*repmat(t,[dev.proj,1]),dt)';