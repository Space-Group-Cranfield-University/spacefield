function R = getRadarMeasurementNoiseCovariance(xTrg, xObs, D_t, SensorParameters)
    dr = norm(xTrg(1:3) - xObs(1:3));
    if SensorParameters.uncertaintySNRdependent
        SNR = computeRadarSNR(1000*dr, D_t, SensorParameters);
        [sigma_r, sigma_v, sigma_theta] = getRadarAccuracy(SNR, SensorParameters);
    else
        [sigma_r, sigma_v, sigma_theta] = getRadarAccuracy(50, SensorParameters);
    end
    sigma_r = sigma_r / 1e3;
    sigma_v = sigma_v / 1e3;
    if SensorParameters.isDoppler
        R = diag([sigma_r^2; sigma_v^2; sigma_theta^2; sigma_theta^2]);
    else
        R = diag([sigma_r^2; sigma_theta^2; sigma_theta^2]);
    end
end