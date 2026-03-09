function [sigma_r, sigma_v, sigma_theta] = getRadarAccuracy(SNR, RadarParameters, CONST)
    % References:
    %   - https://www.mathworks.com/help/radar/ug/measurement-accuracy-bias-and-resolution.html
    if nargin < 4
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 3
        RadarParameters = getStandardRadar("tracking");
    end
    if nargin < 2
        SNR = 1;
    end
    SNR = 10^(SNR / 10);
    theta_b = RadarParameters.theta_b;
    B = RadarParameters.B_r;
    lambda = RadarParameters.lambda;
    f_p = RadarParameters.f_p;
    bias = RadarParameters.bias;
    if nargin < 1
        n_p = 1;
    elseif RadarParameters.radarType == "tracking"
        n_p = f_p * RadarParameters.deltaT;
    elseif RadarParameters.radarType == "fence"
        n_p = RadarParameters.deltaT / RadarParameters.t_s;
    end
    sigma_theta = theta_b / (1.5 * sqrt(n_p * SNR)) + bias(3);
    delta_r = CONST.C / (2 * B);
    delta_v = lambda * f_p / (2 * n_p);
    sigma_r = sqrt(2) * delta_r / sqrt(SNR) + bias(1);
    sigma_v = sqrt(6) * delta_v / (2 * pi * sqrt(SNR)) + bias(2);
end