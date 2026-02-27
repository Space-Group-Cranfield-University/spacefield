function [y, h, R] = getTargetMeasurements(timestep, kTarget, TRG, OBS, ObsTrgCrossVisibilityMat)
    y = [];
    h = @(x) [];
    R = [];
    for kSensor = 1:size(OBS, 2)
        if ObsTrgCrossVisibilityMat(kSensor, kTarget)
            hTemp = @(x) getRaDecFromECI(x);
            R_temp = getOpticalMeasurementNoiseCovariance(OBS(kSensor).SensorParameters);
            y = [y, getMeasurement(TRG(kTarget).xMat(timestep, :)', hTemp, R_temp)];
            h = @(x) [h(x); hTemp(x)];
            R = blkdiag(R, R_temp);
        end
    end
end