function n = sampleNoiseFromCovariance(R)
    n = mvnrnd(zeros(size(R,1), 1), R);
    n = n';
end