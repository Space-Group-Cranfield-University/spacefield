function X = getSigmaPointsUKF(xPrev, P_prev, filterParameters)
    S_prev = chol(P_prev);
    d = size(xPrev, 1);
    lambda = filterParameters.lambda;
    k = sqrt(d+lambda);
    X = [xPrev, xPrev + k*S_prev, xPrev - k*S_prev];
end