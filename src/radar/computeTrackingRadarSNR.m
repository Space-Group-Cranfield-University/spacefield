function SNR = computeTrackingRadarSNR(deltaT, dr, D_t, RadarParameters, CONST)
    % Computes SNR in dB.
    % References:
    %   - Development of the first Portuguese radar tracking sensor for
    %   Space Debris, https://arxiv.org/abs/2102.10457
    if nargin < 5
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 4
        RadarParameters = getStandardRadar("tracking");
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
    sigma = D_t^2;
    k = CONST.K_B;
    n_p = f_p * deltaT;
    P_r = n_p * E_n * P_av * lambda^2 * sigma / ((4*pi)^3 * dr^4);
    N = k * T * B * tau * f_p;
    SNR = P_r / N;
    SNR = 10*log10(SNR) + 2*G + 4*F - N_f - L_s;
end