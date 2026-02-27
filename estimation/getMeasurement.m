function y = getMeasurement(x, h, R)
    n = sampleNoiseFromCovariance(R);
    y = h(x) + n;
end