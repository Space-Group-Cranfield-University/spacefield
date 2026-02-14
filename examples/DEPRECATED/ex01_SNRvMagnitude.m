% Retrieves the dependence of the optical SNR on the apparent magnitude of
% the target object. 
% It uses the formula from Kartunnen et al., "Fundamental Astronomy"
%
% This script replicates results shown in Fig. 2 at https://doi.org/10.1016/j.actaastro.2025.02.019

close all
clear
clc

F_spec = 1/4;
F_diff = @(x) 2/(3*pi)*((pi-x)*cos(x)+sin(x));
rho = 0.3; % Target reflectivity. Are reflectance and geometric albedo the same?
a = 149597870.7; % 1 AU in km
m_sun = -26.73; % visual (i.e. apparent) magnitude Sun (notice: the absolute is 4.83)
r = a; % distance of target from the Sun in km
N_runs = 500;
CONST = initialiseAstronomicalConstants();
%CONST.LAMBDA_VISIBLE = 500*1e-9; % [m], average wavelength visible light
%CONST.C = 299792458;
%CONST.H_SOLAR = 1361; % [W/m2], Solar constant
%CONST.H_PLANCK = 6.62607015*1e-34; % [J/s], Planck's constant
%CONST.NI_VISIBLE = CONST.C /CONST.LAMBDA_VISIBLE; % Average frequency visible light
CAMERA = getStandardCamera(CONST);
CAMERA.D = CAMERA.D_max;
TARGET.rho = rho;
tau = 25*1e-3; % shutter time
minSNR = 10*log10(5);

%% Paper case
[m_target, SNR] = generateMagnitudeSNR(N_runs, m_sun, rho, a, r, tau, F_diff, ...
                F_spec, CAMERA, TARGET, CONST);

figure()
scatter(m_target, SNR, "filled")
hold on
plot([-20 40], [minSNR minSNR], "red")
text(20, 10,"SNR = 5",'Color','red', 'LineWidth', 2, 'FontSize', 14);
fontsize(gca, 12, 'points')
xlabel("apparent magnitude", "FontSize", 14)
ylabel("SNR [dB]", "FontSize", 14)
%title("limiting magnitude", "FontSize", 12)

%% Star tracker case (See Filho et al. 2023)

tau = tau*2;
rho = 0.175;
TARGET.rho = rho;
CAMERA.f = 0.035;
CAMERA.F = 1.6;
CAMERA.F_min = CAMERA.F;
CAMERA.D = CAMERA.f/CAMERA.F;
CAMERA.R_N = 11.97;
CAMERA.I_dark = 1.095;
res = sqrt(1e7); % camera is 10 Mpx
CAMERA.pixel_size = 1.67*1e-6; % Aptina MT9J001 spreadsheet
CAMERA.b(1) = CAMERA.pixel_size * res;
CAMERA.b(2) = CAMERA.b(1);
CAMERA.FOV = deg2rad(10);

[m_target_2, SNR_2] = generateMagnitudeSNR(N_runs, m_sun, rho, a, r, tau, F_diff, F_spec, CAMERA, TARGET, CONST);

figure()
scatter(m_target_2, SNR_2, "filled")
xlabel("apparent magnitude")
ylabel("SNR [dB]")
title("star tracker")

figure()
scatter(m_target, SNR, "filled")
hold on
scatter(m_target_2, SNR_2, "filled")
legend(["camera","star tracker"])
xlabel("apparent magnitude")
ylabel("SNR [dB]")
title("comparison")

%% Functions

function [m_target, SNR] = generateMagnitudeSNR(N_runs, m_sun, rho, a, r, tau, F_diff, F_spec, CAMERA, TARGET, CONST)
    for k = 1:N_runs
        % Generate target parameters
        R = 10.^(-5+rand*3)/2; % Target radius (1 cm to 10 m)
        delta = 10.^(rand*6); % distance of target from observer in km
        TARGET.size = R*1e3; % size in m
        dr = delta*1e3; % relative distance in m
        phi = rand*pi; % phase angle (angle between Sun and Camera wrt target)
        %phi = 0;
        beta = rand;
        %beta = 1;
        v_ta = rand*14*1000; % relative velocity in m/s
        
        % Compute magnitude
        V = m_sun - 2.5*log10(rho*R^2/a^2); % absolute magnitude target
        m_target(k) = V + 5*log10(r*delta/a^2)-2.5*log10(beta*F_diff(phi)+(1-beta)*F_spec); % apparent magnitude target
        
        % Compute SNR (limiting magnitude at SNR = 5)
        SNR(k) = 10*log10(getOpticalSNR(tau, dr, v_ta, CAMERA, TARGET, CONST));
    end
end