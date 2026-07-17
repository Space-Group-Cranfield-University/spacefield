function filterParameters = getStandardParametersUKF()
    %filterParameters.alpha = 1e-3;
    filterParameters.alpha = 0.05;
    filterParameters.beta = 2;
    filterParameters.kappa = 0;
    filterParameters.jitter = 1e-12;
    filterParameters.rThreshold = 10;
    filterParameters.vThreshold = filterParameters.rThreshold / 1000;
end