function CAMERA = getSpaceBasedSSTsensor(sensorName)

    % Function to generate space-based SST sensor parameters.
    % Assumes square FOV.
    % Based on:
    % 
    %   - Ackermann M. et al., 2015, A SYSTEMATIC EXAMINATION OF 
    %   GROUND-BASED AND SPACE-BASED APPROACHES TO OPTICAL DETECTION AND 
    %   TRACKING OF SATELLITES
    %   Available from: 
    %   https://www.spacesymposium.org/wp-content/uploads/2017/10/M.Ackermann_31st_Space_Symposium_Tech_Track_paper.pdf
   
    if nargin < 1
        sensorName = "";
    end

    switch sensorName
        case "Sapphire"
            CAMERA = getSpaceBasedSST_Sapphire();
        otherwise
            CAMERA = getSpaceBasedSST_SBSS();
    end

end

function CAMERA = getSpaceBasedSST_Sapphire()
    lambda = initializeAstronomicalConstants().LAMBDA_VISIBLE; % half-band wavelength
    CAMERA.D = 0.15;                  % Aperture
    CAMERA.pixel_size = 12*1e-6;        % Pixel size
    CAMERA.pixel_resolution = sqrt(2*1e6); % Resolution in pixels
    CAMERA.FOV = deg2rad(1.4); % FIeld of view
    CAMERA.Q = 0.6;                     % Quantum efficiency
    %CAMERA.f = 0.0518;                  % Focal length
    CAMERA.R_N_sq = 2;                 % Readout noise
    CAMERA.I_dark = 0.1;               % Dark current
    CAMERA.C = 8*1e4;                   % Pixel well depth
    CAMERA.K = 0.1;                     % Sky count rate
    CAMERA.tau = 0.1;                   % Exposure time
    CAMERA.saturation_threshold = 0.9;  % Saturation threshold
    CAMERA.SNR_threshold = 5;           % SNR threshold for positive detection
    CAMERA.std_attitude = 8 / 206265;   % Inaccuracy due to attitude estimation

    % Parameter processing
    CAMERA.b = CAMERA.pixel_resolution*CAMERA.pixel_size;
    CAMERA.f = CAMERA.b / (2 * tan(CAMERA.FOV / 2) );
    %CAMERA.FOV = 2*atan(CAMERA.b / (2*CAMERA.f));
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

function CAMERA = getSpaceBasedSST_SBSS()
    
    lambda = initializeAstronomicalConstants().LAMBDA_VISIBLE; % half-band wavelength
    CAMERA.D = 0.3;                  % Aperture
    CAMERA.pixel_size = 12*1e-6;        % Pixel size
    CAMERA.pixel_resolution = sqrt(2*1e6); % Resolution in pixels
    CAMERA.FOV = deg2rad(sqrt(8)); % FIeld of view
    CAMERA.Q = 0.6;                     % Quantum efficiency
    %CAMERA.f = 0.0518;                  % Focal length
    CAMERA.R_N_sq = 2;                 % Readout noise
    CAMERA.I_dark = 0.1;               % Dark current
    CAMERA.C = 8*1e4;                   % Pixel well depth
    CAMERA.K = 0.1;                     % Sky count rate
    CAMERA.tau = 0.1;                   % Exposure time
    CAMERA.saturation_threshold = 0.9;  % Saturation threshold
    CAMERA.SNR_threshold = 5;           % SNR threshold for positive detection
    CAMERA.std_attitude = 8 / 206265;   % Inaccuracy due to attitude estimation

    % Parameter processing
    CAMERA.b = CAMERA.pixel_resolution*CAMERA.pixel_size;
    CAMERA.f = CAMERA.b / (2 * tan(CAMERA.FOV / 2) );
    %CAMERA.FOV = 2*atan(CAMERA.b / (2*CAMERA.f));
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