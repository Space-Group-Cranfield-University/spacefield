function R = getOpticalMeasurementNoiseCovariance(SensorParameters)
    R = eye(2) * SensorParameters.sigma^2;
end