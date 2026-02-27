function OBS = setConstellationSensors(OBS, SensorParameters)
    if nargin < 1
        SensorParameters = getReducedSensorParameters;
    end
    for k = 1:size(OBS, 2)
        OBS(k).SensorParameters = SensorParameters;
    end
end