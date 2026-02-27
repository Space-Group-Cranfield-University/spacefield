function [xNext, P_next] = stepUKF(deltaT, xPrev, P_prev, y, f, h, isScaled, D)
    if nargin < 7
        isScaled = 1;
        D = getStandardScaleMatrix;
    end
    if ~isScaled
        D = eye(size(xPrev, 2));
    end
    [xPrev, P_prev] = scaleInputs(xPrev, P_prev);
    [wMeanVec, wCovVec, filterParameters] = getWeightsUKF(size(xPrev, 1));
    [xPred, P_pred, X_pred, X_var] = ...
        predictUT(deltaT, xPrev, P_prev, f, wMeanVec, wCovVec, filterParameters, isScaled, D, 1);
    Y_pred = [];
    for k = 1:size(X_pred, 2)
        Y_pred = [Y_pred, h(D \ X_pred(:, k))];
    end
    yPred = Y_pred * wMeanVec;
    Y_var = Y_pred - yPred;
    Pyy = Y_var * diag(wCovVec) * Y_var';
    Pyy = 0.5*(Pyy + Pyy');   % symmetry fix
    Pyy = Pyy + filterParameters.jitter*eye(size(Pyy, 1));
    Pxy = X_var * diag(wCovVec) * Y_var';
    L = chol(Pyy,'lower'); % Helps numerical robustness
    K = (Pxy / L') / L;
    xNext = xPred + K * (y - yPred);
    P_next = P_pred - K * Pyy * K';
    P_next = ( P_next + P_next' ) / 2; % Helps positive definitness
    [xNext, P_next] = unscaleOutputs(xNext, P_next);
end