function [xPred, P_pred, X_pred, X_var] = ...
        predictUT(deltaT, xPrev, P_prev, f, wMeanVec, wCovVec, filterParameters, isScaled, D, calledWithinFilter)
    if nargin < 10
        calledWithinFilter = 0;
    end
    if nargin < 8
        isScaled = 1;
        D = getStandardScaleMatrix;
    end
    if nargin < 5
        [wMeanVec, wCovVec, filterParameters] = getWeightsUKF(size(xPrev, 1));
    end
    if ~isScaled
        D = eye(size(xPrev, 2));
    end
    if ~calledWithinFilter
        [xPrev, P_prev] = scaleInputs(xPrev, P_prev, D);
    end
    X = getSigmaPointsUKF(xPrev, P_prev, filterParameters);
    X_pred = [];
    for k = 1:size(X, 2)
        X_pred = [X_pred, D*f(D \ X(:, k))];
    end
    xPred = X_pred * wMeanVec;
    X_var = X_pred - xPred;
    Q = getProcessNoiseCovariance(deltaT);
    P_pred = X_var * diag(wCovVec) * X_var' + Q;
    if ~calledWithinFilter
        [xPred, P_pred] = unscaleOutputs(xPred, P_pred, D);
    end
    P_pred = (P_pred + P_pred') / 2; % Helps positive definitness
end