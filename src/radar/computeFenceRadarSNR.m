function SNR = computeFenceRadarSNR(deltaT, dr, D_t, RadarParameters, CONST)
    if nargin < 5
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 4
        RadarParameters = getStandardRadar("fence");
    end
    P_av = RadarParameters.P_av;
    f_p = RadarParameters.f_p;
    G = RadarParameters.G; % dB
    lambda = RadarParameters.lambda;
    E_n = RadarParameters.E_n;
    F = RadarParameters.F; % dB
    tau = RadarParameters.tau;
    N_f = RadarParameters.N_f; % dB
    T = RadarParameters.T;
    B = RadarParameters.B_r;
    L_s = RadarParameters.L_s; % dB
    Om_s = RadarParameters.Om_s;
    sigma = D_t^2;
    k = CONST.K_B;
    % I use deltaT instead of t_s as I do multiple scans over a given
    % interval
    P_r = E_n * P_av * lambda^2 * sigma * deltaT / (16*pi^2 * dr^4 * Om_s);
    N = k * T * B * tau * f_p;
    SNR = P_r / N;
    SNR = 10*log10(SNR) + G + 4*F - N_f - L_s;

end