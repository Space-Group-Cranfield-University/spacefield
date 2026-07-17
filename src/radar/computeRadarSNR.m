function SNR = computeRadarSNR(dr, D_t, RadarParameters, CONST)
    % Computes SNR in dB
    if nargin < 5
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 4
        RadarParameters = getStandardRadar;
    end
    if RadarParameters.radarType == "tracking"
        SNR = computeTrackingRadarSNR(RadarParameters.deltaT, dr, D_t, RadarParameters, CONST);
    elseif RadarParameters.radarType == "fence"
        SNR = computeFenceRadarSNR(RadarParameters.deltaT, dr, D_t, RadarParameters, CONST);
    end
    if SNR > RadarParameters.maxSNR
        SNR = RadarParameters.maxSNR;
    end
end