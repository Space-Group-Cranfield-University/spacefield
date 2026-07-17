function TRG = updateAllTargetStates(timestep, deltaT, TRG, OBS, ObsTrgCrossVisibilityMat)
    for kTarget = 1:size(TRG, 2)
        [y, h, R] = getTargetMeasurements(timestep, kTarget, TRG, OBS, ObsTrgCrossVisibilityMat);
        TRG = estimateTargetState(deltaT, timestep, kTarget, TRG, y, h);
        TRG = checkThresholds(timestep, kTarget, TRG);
    end
end