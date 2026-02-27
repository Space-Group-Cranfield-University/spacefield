function [xNext, P_next] = stepUKF(xPrev, P_prev, y, f, h)
    [wMeanVec, wCovVec, filterParameters] = getWeightsUKF(size(xPrev, 1));
    [xPred, P_pred, X_pred, X_var] = ...
        predictUT(xPrev, P_prev, f, wMeanVec, wCovVec, filterParameters);
    Y_pred = [];
    for k = 1:size(X_pred, 2)
        Y_pred = [Y_pred, h(X_pred(:, 2))];
    end
    yPred = Y_pred * wMeanVec;
    Y_var = Y_pred - yPred;
    Pyy = Y_var * diag(wCovVec) * Y_var';
    Pxy = X_var * diag(wCovVec) * Y_var';
    K = Pxy * inv(Pyy);
    xNext = xPred + K * (y - yPred);
    P_next = P_pred - K * Pyy * K';
end