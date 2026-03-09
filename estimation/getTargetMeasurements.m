function [y, h, R] = getTargetMeasurements(timestep, kTarget, TRG, OBS, ObsTrgCrossVisibilityMat)
    % Deprecated (yielded singular covariance matrices)
    y = [];
    h = @(x) [];
    R = [];
    for kSensor = 1:size(OBS, 2)
        if ObsTrgCrossVisibilityMat(kSensor, kTarget)
            if OBS(kSensor).SensorParameters.sensorType == "optical"
                hTemp = @(x) getOpticalMeasurement(x, OBS(kSensor).xMat(timestep, :)');
                R_temp = getOpticalMeasurementNoiseCovariance(OBS(kSensor).SensorParameters);
            elseif OBS(kSensor).SensorParameters.sensorType == "radar"
                hTemp = @(x) getRadarMeasurement...
                    (x, OBS(kSensor).xMat(timestep, :)', OBS(kSensor).SensorParameters);
                R_temp = ...
                    getRadarMeasurementNoiseCovariance...
                    (TRG(kTarget).xMat(timestep, :)', ...
                    OBS(kSensor).xMat(timestep, :)', ...
                    TRG(kTarget).D_t, OBS(kSensor).SensorParameters);
            end
            y = [y; getMeasurement(TRG(kTarget).xMat(timestep, :)', hTemp, R_temp)];
            h = @(x) [h(x); hTemp(x)];
            R = blkdiag(R, R_temp);
        end
    end
end

function [y, h, R] = getTargetMeasurementsSingle(timestep, kTarget, TRG, OBS, ObsTrgCrossVisibilityMat)
    R = getOpticalMeasurementNoiseCovariance(OBS(1).SensorParameters);
    N_sensors = sum(ObsTrgCrossVisibilityMat(:, kTarget));
    if N_sensors > 0
        R = R / N_sensors;
        h = @(x) getRaDecFromECI(x);
        y = getMeasurement(TRG(kTarget).xMat(timestep, :)', h, R);
    else
        y = [];
        h = @(x) [];
        R = [];
    end
end