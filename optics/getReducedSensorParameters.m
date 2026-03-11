function SensorParameters = getReducedSensorParameters(FOV, sensorType, radarType)
    if nargin < 2
        sensorType = "optical";
    end
    if nargin < 1
        FOV = deg2rad(8.1);
    end
    if sensorType == "radar" && nargin == 3
        SensorParameters = getStandardRadar(radarType);
    elseif sensorType == "radar"
        warning("Radar requires additional argument")
    else
        SensorParameters.V_lim = 14.6;
        SensorParameters.alpha_e = deg2rad(15);
        SensorParameters.alpha_s = SensorParameters.alpha_e;
        SensorParameters.FOV = FOV;
        SensorParameters.halfFOV = SensorParameters.FOV / 2;
        SensorParameters.sigma = convertArcsecToRad(8.3);
        SensorParameters.sensorType = "optical";
        SensorParameters.R_H = 6471;
        SensorParameters.maxSlewAngle = getStandardGreedyOptions().maxSlewAngle;
        SensorParameters.minEl = getStandardGreedyOptions().minEl;
        SensorParameters.maxEl = getStandardGreedyOptions().maxEl;
        SensorParameters.orbitalFractionWellTracked = 0.5;
    end
end