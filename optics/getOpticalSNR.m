function SNR = getOpticalSNR(tau, dr_obj, v_ta, CAMERA, TRG, CONST)

    % Function to get measurement parameters such as accuracy
    % Camera coordinate system is centred at the optical focus and with z
    % perpendicular to image plane.
    % Based on:
    % Pineau D & Felicetti L (2023) Design of an optical system for a 
    % Multi-CubeSats debris surveillance mission, Acta Astronautica, 
    % 210 (September) 535-546.
    %
    % Assumption: measurements are taken at tau/2

    Q = CAMERA.Q; % Quantum efficiency
    H = CONST.H_SOLAR; % Solar constant
    h = CONST.H_PLANCK; % Planck constant
    ni = CONST.NI_VISIBLE; % Average frequency visible light
    rho = TRG.rho; % Target reflectivity
    r = TRG.size; % Target size
    D = CAMERA.D; % Camera aperture diameter
    % tau = CAMERA.tau; % Exposure time of image, set by user
    d = norm(dr_obj);
    lambda = CONST.LAMBDA_VISIBLE; % Average wavelength visible light
    f = CAMERA.f; % [m], focal length
    p = CAMERA.pixel_size; % [m], pixel size
    K = CAMERA.K; % sky count rate
    I_dark = CAMERA.I_dark; % dark current
    R_N = CAMERA.R_N; % readout noise
    
    gamma = 0.4; % average value, does NOT consider object attitude, 
    % it is the same results as if considering a spherical target 
    % and average camera positioning

    S = Q * ( H / ( 8 * h * ni) ) * ( rho * r^2 * gamma ) * ( D^2 * tau / d^2 );
    
    n_ip = ( 2.44 * lambda ) * ( v_ta * tau * f^2 ) / ( dr_obj(1) * D * p^2 ); % number of illuminated pixels
    N = sqrt( S + n_ip * ( ( K + I_dark ) * tau + R_N^2 ) );

    SNR = S / N; % NOTICE: for multiple exposures SNR ~ SNR*sqrt(n_exposures) https://www.stsci.edu/instruments/wfpc2/Wfpc2_hand_current/ch6_exposuretime6.html

end