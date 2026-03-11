function RadarParameters = getStandardRadar(radarType)
    % References:
    %   - Development of the first Portuguese radar tracking sensor for
    %   Space Debris, https://arxiv.org/abs/2102.10457
    %   - Radial velocity ambiguity resolution for the pulsed space
    %   surveillance radar GESTRA, 
    %   https://conference.sdo.esoc.esa.int/proceedings/sdc9/paper/236
    %   - https://en.wikipedia.org/wiki/AN/FPS-16_Instrumentation_Radar
    if nargin < 1
        radarType = "tracking";
    end
    CONST = initializeAstronomicalConstants;
    RadarParameters.sensorType = "radar";
    RadarParameters.radarType = radarType;
    RadarParameters.G = 60; % Gain [dB]
    RadarParameters.theta_b = deg2rad(0.2); % beamwidth [rad], 12 arcmin
    % Notice: local reference systems has z axis pointing down
    RadarParameters.maxEl = -deg2rad(20); % Minimum elevation [rad]
    RadarParameters.minEl = -deg2rad(90); % Minimum elevation [rad]
    RadarParameters.maxElRate = deg2rad(2); % Elevation rate [rad/s]
    RadarParameters.G_T = 40; % Gain over temperature G/T [dB]
    RadarParameters.L_s = 1; % System losses [dB]
    RadarParameters.E_n = 1; % Pulse integration efficiency [-]
    RadarParameters.F = 0; % Propagation effects [dB]
    RadarParameters.N_f = 1; % Receiver noise figure [dB]
    RadarParameters.T = 288; % Receiver temperture [K]
    if radarType == "tracking"
        RadarParameters.P_t = 5*1e5; % Max transmitter power [W]
        RadarParameters.f = 1e10; % Nominal frequency [Hz], set to X-band
        RadarParameters.maxAzRate = deg2rad(10); % Azimuth rate [rad/s]
        RadarParameters.f_p = 300; % Pulse repetition frequency [Hz]
        RadarParameters.tau = 1 * 1e-3; % Pulse width [s]
        RadarParameters.D = 10; % Aperture [m], this should not be used as only the effective area goes into the radar equation
        RadarParameters.B_r = 1e6; % Receiver bandwidth [Hz]
    elseif radarType == "fence"
        RadarParameters.P_t = 5*1e5; % Max transmitter power [W]
        RadarParameters.f = 5*1e8; % Nominal frequency [Hz], set to UHF
        RadarParameters.azRange = deg2rad(60); % Fence width
        RadarParameters.elRange = deg2rad(2); % Fence width
        RadarParameters.f_p = 150; % Pulse repetition frequency [Hz]
        RadarParameters.tau = 5 * 1e-3; % Pulse width [s]
        RadarParameters.B_r = 1e5; % Receiver bandwidth [Hz]
        RadarParameters.Om_s = ...
            pi * RadarParameters.azRange * RadarParameters.elRange; % Scanned solid angle [srad]
        RadarParameters.t_s = ...
            RadarParameters.Om_s / ...
            (RadarParameters.theta_b * RadarParameters.f_p); % Scanning time for a single pass, not used in radar equation
    end
    RadarParameters.P_av = RadarParameters.P_t * RadarParameters.tau * RadarParameters.f_p;
    RadarParameters.lambda = CONST.C / RadarParameters.f;
    RadarParameters.bias = [10; 0.1; 10/206000];
    RadarParameters.maxSNR = 50; % [dB]
    RadarParameters.deltaT = 1; % Integration time [s]
    RadarParameters.minSNR = 5; % [dB]
    RadarParameters.FOV = deg2rad(40);
    RadarParameters.halfFOV = RadarParameters.FOV / 2;
    RadarParameters.R_H = 6371;
    RadarParameters.alpha_e = -RadarParameters.maxEl;
    RadarParameters.maxSlewAngle = deg2rad(10);
    RadarParameters.isDoppler = 0;
    RadarParameters.uncertaintySNRdependent = 1;
    RadarParameters.orbitalFractionWellTracked = 1;
end