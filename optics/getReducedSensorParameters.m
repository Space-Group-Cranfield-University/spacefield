function SensorParameters = getReducedSensorParameters()
    SensorParameters.V_lim = 14.6;
    SensorParameters.alpha_e = deg2rad(25);
    SensorParameters.alpha_s = SensorParameters.alpha_e;
    SensorParameters.FOV = deg2rad(8.1);
    SensorParameters.halfFOV = SensorParameters.FOV / 2;
end