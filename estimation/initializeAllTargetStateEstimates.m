function TRG = initializeAllTargetStateEstimates(TRG, timeVec, P_0)
    sigma_r = 1;
    sigma_v = sigma_r / 1000;
    if nargin < 3
        P_0 = diag( ...
            [sigma_r^2, sigma_r^2, sigma_r^2, ...
            sigma_v^2, sigma_v^2, sigma_v^2]);
    end
    for k = 1:size(TRG, 2)
        TRG(k).xMatEst = zeros(size(timeVec, 2), 6);
        TRG(k).P = zeros(6, 6, size(timeVec, 2));
        dx = sampleNoiseFromCovariance(P_0);
        TRG(k).xMatEst(1, :) = TRG(k).x0' + dx';
        TRG(k).P(:,:,1) = P_0;
    end
end