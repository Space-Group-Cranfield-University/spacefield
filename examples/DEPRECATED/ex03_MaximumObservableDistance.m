% This script computes the maximum observable distance as a function of the
% target size.
%
% This script replicates results shown in Fig. 6 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

tau = 0.025;
dr_obj = 100*1000;
v_ta = 10.2*1000; % [m/s]
CONST = initialiseAstronomicalConstants();
CONST.LAMBDA_VISIBLE = 500*1e-9; % [m], average wavelength visible light
CONST.C = 299792458;
CONST.H_SOLAR = 1361; % [W/m2], Solar constant
CONST.H_PLANCK = 6.62607015*1e-34; % [J/s], Planck's constant
CONST.NI_VISIBLE = CONST.C /CONST.LAMBDA_VISIBLE; % Average frequency visible light
CAMERA = getStandardCamera(CONST);
CAMERA.D = CAMERA.D_max;
TRG.rho = 0.3;
size = 1e-3; % [m]
SNR_threshold = 10*log10(5); % [dB]
FOV = CAMERA.FOV;
tanFOV2 = tan(FOV/2);
h = 400; % Constellation altitude
R = CONST.R_E + h;

%% Find SNR = 7.5 dB curve in (size, dr) space

N = 100;
size_vect = 10.^linspace(-3,0,N);
dr_0 = 1e4; % [m]

for k = 1:N
    size_temp = size_vect(k);
    dr_vect(k) = fsolve(@(dr) getSNRdB(tau, dr, v_ta, size_temp, CAMERA, TRG, CONST) - SNR_threshold, dr_0);
end
dr_min = 2*v_ta*tau/tanFOV2;
difference = dr_vect/1000 - dr_min * ones(1,N) / 1000;
indexVec = find(difference > 0);

figure()
loglog(size_vect*1000, dr_vect/1000, "LineWidth", 2)
hold on
loglog(size_vect*1000, dr_min * ones(1,N) / 1000, "r", 'LineWidth', 2)
x = [size_vect(indexVec)*1000, fliplr(size_vect(indexVec)*1000)];
dr_min_vec = dr_min * ones(1,N) / 1000;
inBetween = [dr_min_vec(indexVec), fliplr(dr_vect(indexVec)/1000)];
flipped = fliplr(inBetween);
fill(x, flipped, 'y');
text(size_vect(ceil(end/3-1))*1000, 500 + dr_vect(end/2)/1000,"SNR = 5",'Color','blue', 'LineWidth', 2, 'FontSize', 14);
text(size_vect(ceil(end/3))*1000, dr_min / 1000 - 3,"transit time = "+char([0xD835 0xDF0F]),'Color','red', 'LineWidth', 2, 'FontSize', 14);
fontsize(gca, 12, 'points')
xlabel("target size [mm]", "FontSize", 14)
ylabel("distance [km]", "FontSize", 14)

%% Functions

function SNR_dB = getSNRdB(tau, dr_obj, v_ta, size_in, CAMERA, TRG, CONST)

    Q = CAMERA.Q; % Quantum efficiency
    H = CONST.H_SOLAR; % Solar constant
    h = CONST.H_PLANCK; % Planck constant
    ni = CONST.NI_VISIBLE; % Average frequency visible light
    rho = TRG.rho; % Target reflectivity
    r = size_in; % Target size
    D = CAMERA.D; % Camera aperture diameter
    % tau = CAMERA.tau; % Exposure time of image, set by user
    d = dr_obj;
    lambda = CONST.LAMBDA_VISIBLE; % Average wavelength visible light
    f = CAMERA.f; % [m], focal length
    p = CAMERA.pixel_size; % [m], pixel size
    K = CAMERA.K; % sky count rate
    I_dark = CAMERA.I_dark; % dark current
    R_N = CAMERA.R_N; % readout noise
    
    gamma = 0.4; %(average value, does NOT consider object attitude, it is the same results as if considering a spherical target and average camera positioning. ISSUE: camera gets irradiated by LESS light than that captured by the target, which is 2*P/pi ~ 0.64*P!!!)
    % gamma = cos(theta_s) * cos(theta_c);
    S = Q * ( H / ( 8 * h * ni) ) * ( rho * r^2 * gamma ) * ( D^2 * tau / d^2 );
    
    n_ip = ( 2.44 * lambda ) * ( v_ta * tau * f^2 ) / ( dr_obj * D * p^2 ); % number of illuminated pixels
    N = sqrt( S + n_ip * ( ( K + I_dark ) * tau + R_N^2 ) );

    SNR = S / N; % NOTICE: for multiple exposures SNR ~ SNR*sqrt(n_exposures) https://www.stsci.edu/instruments/wfpc2/Wfpc2_hand_current/ch6_exposuretime6.html

    SNR_dB = 10*log10(SNR);
end