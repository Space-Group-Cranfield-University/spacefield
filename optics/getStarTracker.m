function CAMERA = getStarTracker(starTracker)

    % Function to generate Star Tracker sensor parameters.
    % Assumes square FOV.
    % Based on AD-1 Star Tracker (Mars Bureau) from:
    % 
    %   - SD. On-Orbit Resident Space Object (RSO) Detection Using 
    %   Commercial Grade Star Trackers; 2019. 
    %   Available from: http://hdl.handle.net/10315/36799
    %   - Filho et al., 2023, Satellite Star Tracker Breadboard with Space 
    %   Debris Detection Capability for LEO
    %   Available from: 
    %   https://iopscience.iop.org/article/10.1088/1742-6596/2526/1/012119
    
    if nargin < 1
        starTracker = "";
    end

    switch starTracker
        case "Filho"
            CAMERA = getStarTrackerFilho();
        otherwise
            CAMERA = getStarTrackerClemens();
    end

end

function CAMERA = getStarTrackerFilho()
    CAMERA.F = 1.6;
    CAMERA.FOV = deg2rad( 10.5 + 7.5 ) / 2;
    CAMERA.D = 0.022;
    CAMERA.R_N_sq = 11.97;
    CAMERA.I_dark = 1.095;
    CAMERA.tau = 0.1;
    lambda = 550*1e-9;                  % half-band wavelength
    CAMERA.Q = 0.2;
    CAMERA.C = 8*1e4;                   % Pixel well depth
    CAMERA.K = 3.7;                     % Sky count rate
    CAMERA.saturation_threshold = 0.9;  % Saturation threshold
    CAMERA.SNR_threshold = 5;           % SNR threshold for positive detection
    CAMERA.std_attitude = 8 / 206265;   % Inaccuracy due to attitude estimation

    % Data processing
    CAMERA.f = CAMERA.D * CAMERA.F;
    CAMERA.b = 2 * tan(CAMERA.FOV / 2) * CAMERA.f;
    CAMERA.pixel_resolution = 512;
    CAMERA.pixel_size = CAMERA.b / CAMERA.pixel_resolution;
    CAMERA.halfFOV = CAMERA.FOV / 2;
    CAMERA.FOV_eq = getSquareToConicalFOV(CAMERA.FOV);
    CAMERA.halfFOV_eq = CAMERA.FOV_eq / 2;
    
    % Compute angular accuracies (the higher the accuracy, the lower the
    % number). Given in [rad]
    CAMERA.std_pixel = CAMERA.pixel_size / ( sqrt(12) * CAMERA.f );
    CAMERA.std_airy = 1.22 * lambda / CAMERA.D;
    CAMERA.std = sqrt( CAMERA.std_pixel^2 + CAMERA.std_airy^2 + CAMERA.std_attitude^2 );
end

function CAMERA = getStarTrackerClemens()
    lambda = 700*1e-9;                  % half-band wavelength
    CAMERA.D = 0.0288;                  % Aperture
    CAMERA.pixel_size = 23*1e-6;        % Pixel size
    CAMERA.pixel_resolution = 512;      % Resolution in pixels
    CAMERA.Q = 0.2;                     % Quantum efficiency
    CAMERA.f = 0.0518;                  % Focal length
    CAMERA.R_N_sq = 40;                 % Readout noise
    CAMERA.I_dark = 3000;               % Dark current
    CAMERA.C = 8*1e4;                   % Pixel well depth
    CAMERA.K = 3.7;                     % Sky count rate
    CAMERA.tau = 0.5;                   % Exposure time
    CAMERA.saturation_threshold = 0.9;  % Saturation threshold
    CAMERA.SNR_threshold = 5;           % SNR threshold for positive detection
    CAMERA.std_attitude = 8 / 206265;   % Inaccuracy due to attitude estimation

    % Parameter processing
    CAMERA.b = CAMERA.pixel_resolution*CAMERA.pixel_size;
    CAMERA.FOV = 2*atan(CAMERA.b / (2*CAMERA.f));
    CAMERA.halfFOV = CAMERA.FOV / 2;
    CAMERA.FOV_eq = getSquareToConicalFOV(CAMERA.FOV);
    CAMERA.halfFOV_eq = CAMERA.FOV_eq / 2;
    CAMERA.F = CAMERA.f / CAMERA.D;
    
    % Compute angular accuracies (the higher the accuracy, the lower the
    % number). Given in [rad]
    CAMERA.std_pixel = CAMERA.pixel_size / ( sqrt(12) * CAMERA.f );
    CAMERA.std_airy = 1.22 * lambda / CAMERA.D;
    CAMERA.std = sqrt( CAMERA.std_pixel^2 + CAMERA.std_airy^2 + CAMERA.std_attitude^2 );
end