function CAMERA = getStandardCamera(CONST)

    % Gives optical sensor parameters for a standard case. The standard
    % case is an off-the-shelf (realitvely cheap) camera with realistic 
    % parameters such as exposure time, aperture, focal length.
    %
    % Based on:
    %   Pineau D & Felicetti L (2023) Design of an optical system for a 
    %   Multi-CubeSats debris surveillance mission, Acta Astronautica, 
    %   210 (September) 535-546.
    %
    % The paper considers telescopic cameras such as Canon EF 400mm f/2.8L IS
    % II USM, whose price is around USD 16000.
    %
    % NOTICE: FOV does NOT depend on aperture!!!

    if nargin < 1
        lambda = initialiseAstronomicalConstants().LAMBDA_VISIBLE;
    else
        lambda = CONST.LAMBDA_VISIBLE;
    end
    
    CAMERA.b = [0.05, 0.05]; % [m], sensor sizes b_1, b_2
    CAMERA.pixel_size = 12*1e-6; % [m], pixel size
    CAMERA.Q = 0.61; % quantum efficiency
    CAMERA.C = 8*1e4; % pixel well depth in number of electrons
    CAMERA.R_N_sq = 4.3; % readout noise in number of electrons per pixel
    CAMERA.I_dark = 35; % dark current in number of electrons per second per pixel
    CAMERA.K = 3.7; % sky count rate in number of electrons per second per pixel
    CAMERA.F = 2.8;
    CAMERA.F_min = 2.8; % min focal ratio
    CAMERA.F_max = 32; % max focal ratio
    CAMERA.f = 0.4; % [m], focal length
    CAMERA.tau = 0.01;
    CAMERA.tau_min = 1/120; % [s], min exposure time
    CAMERA.tau_max = 120; % [s], max exposure time
    CAMERA.SNR_threshold = 5; % minimum signal to noise ratio for detectabel images. It is NOT in dB!!!
    CAMERA.saturation_threshold = 0.9; % camera saturates at S + N = n_ip * sat_threshold * C
    CAMERA.std_attitude = 8 / 206265; % Angular resolution due to attitude errors (assuming optical sensor is mounted on a CubeSat system)

    % Post process
    CAMERA.D = CAMERA.f / CAMERA.F;
    CAMERA.D_min = CAMERA.f / CAMERA.F_max; % [m], min aperture
    CAMERA.D_max = CAMERA.f / CAMERA.F_min; % [m], max aperture
    CAMERA.pixels = floor(CAMERA.b / CAMERA.pixel_size);
    CAMERA.FOV = 2*atan(CAMERA.b / (2*CAMERA.f)); % [rad], assuming square field of view. This is not usually the case, the y-axis will have a smaller angle. Notice this is 2*atan(b/(2f))
    CAMERA.pixel_resolution = CAMERA.b(1) / CAMERA.pixel_size;
    
    % Compute angular accuracies (the higher the accuracy, the lower the
    % number). Given in [rad]
    % CAMERA.std_pixel = 2 * CAMERA.pixel_size / ( sqrt(2) * CAMERA.f );
    CAMERA.std_pixel = CAMERA.pixel_size / ( sqrt(12) * CAMERA.f );
    CAMERA.std_airy = 1.22 * lambda / CAMERA.D;
    CAMERA.std_airy_min = 1.22 * lambda / CAMERA.D_max;
    CAMERA.std_airy_max = 1.22 * lambda / CAMERA.D_min;
    %CAMERA.std = max( CAMERA.std_pixel, CAMERA.std_airy );
    CAMERA.std_min = max( CAMERA.std_pixel, CAMERA.std_airy_min );
    CAMERA.std_max = max( CAMERA.std_pixel, CAMERA.std_airy_max );
    CAMERA.std = sqrt( CAMERA.std_pixel^2 + CAMERA.std_airy^2 + CAMERA.std_attitude^2 );
    
end