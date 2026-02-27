function [wMeanVec, wCovVec, filterParameters] = getWeightsUKF(d, filterParameters)
    % Notice that for small alpha, due to numerical error, the sum of the
    % weights is not exactly one.
    if nargin < 2
        filterParameters = getStandardParametersUKF;
    end
    if nargin < 1
        d = 6;
    end
    filterParameters.lambda     =  - d + filterParameters.alpha^2 * ...
                                (d + filterParameters.kappa);
    wMeanVec                    = zeros(2*d + 1, 1);
    wCovVec                     = zeros(2*d + 1, 1);
    wMeanVec(1)                 = filterParameters.lambda / ...
                                (d + filterParameters.lambda);
    wCovVec(1)                  = wMeanVec(1) + ...
                                (1 - filterParameters.alpha^2 ...
                                + filterParameters.beta);
    wMeanVec(2:end)             = 1 / (2* (d + filterParameters.lambda) );
    wCovVec(2:end)              = wMeanVec(2:end);
end