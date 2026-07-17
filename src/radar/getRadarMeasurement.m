function y = getRadarMeasurement(xTrg, xObs, SensorParameters)
    rTrg = xTrg(1:3);
    rObs = xObs(1:3);
    vTrg = xTrg(4:6);
    vObs = xObs(4:6);
    dr = rTrg - rObs;
    dv = vTrg - vObs;
    AzEl = getOpticalMeasurement(xTrg, xObs);
    rho = norm(dr);
    rho_dot = dr' * dv / rho;
    if SensorParameters.isDoppler
        y = [rho; rho_dot; AzEl];
    else
        y = [rho; AzEl];
    end
end