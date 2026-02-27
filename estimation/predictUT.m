function [xPred, P_pred, X_pred, X_var] = predictUT(xPrev, P_prev, f, wMeanVec, wCovVec, filterParameters)
    if nargin < 4
        [wMeanVec, wCovVec, filterParameters] = getWeightsUKF(size(xPrev, 1));
    end
    X = getSigmaPointsUKF(xPrev, P_prev, filterParameters);
    X_pred = [];
    for k = 1:size(X, 2)
        X_pred = [X_pred, f(X(:, k))];
    end
    xPred = X_pred * wMeanVec;
    X_var = X_pred - xPred;
    P_pred = X_var * diag(wCovVec) * X_var';
end