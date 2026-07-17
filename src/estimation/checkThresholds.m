function TRG = checkThresholds(timestep, kTarget, TRG, filterParameters)
    if nargin < 4
        filterParameters = getStandardParametersUKF;
    end
    EE = TRG(kTarget).xMatEst(timestep, :) - TRG(kTarget).xMat(timestep, :);
    for k = 1:3
        if abs(EE(k)) > filterParameters.rThreshold
            EE(k) = sign(EE(k)) * filterParameters.rThreshold;
        end
    end
    for k = 4:6
        if abs(EE(k)) > filterParameters.vThreshold
            EE(k) = sign(EE(k)) * filterParameters.vThreshold;
        end
    end
    TRG(kTarget).xMatEst(timestep, :) = TRG(kTarget).xMat(timestep, :) + EE;
end