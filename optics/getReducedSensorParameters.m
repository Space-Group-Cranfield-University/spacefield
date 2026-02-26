function SensorParameters = getReducedSensorParameters(FOV)
    if nargin < 1
        FOV = deg2rad(8.1);
    end
    SensorParameters.V_lim = 14.6;
    SensorParameters.alpha_e = deg2rad(15);
    SensorParameters.alpha_s = SensorParameters.alpha_e;
    SensorParameters.FOV = FOV;
    SensorParameters.halfFOV = SensorParameters.FOV / 2;
end