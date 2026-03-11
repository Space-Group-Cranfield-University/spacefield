function [nObservedTargets, nWeightedTargets, taggedTargets, targetWeights] = ...
                countTaskingVisibilityOverTargetPopulation...
                (timestep, OBS, TRG, prevObsTrgCrossVisibilityVec, dirSun, ...
                OPTIONS)
    nObservedTargets = 0;
    nWeightedTargets = 0;
    taggedTargets = [];
    targetWeights = [];
    for k = 1:size(TRG, 2)
        % Pointing direction must be converted to ECI to work within
        % visibility function.
        dirPointing = getPointingDirectionInECI(OBS, timestep);
        rTrg = TRG(k).xMat(timestep, 1:3)';
        if OPTIONS.isTargetStateUncertain ...
                && ~(OBS.SensorParameters.sensorType == "radar")
            rTrg = rTrg + randn(3, 1) * OPTIONS.positionError;
        end
        if OPTIONS.isTargetStateUncertain...
                && norm(rTrg) < OPTIONS.positionError * OBS.SensorParameters.FOV ...
                && ~(OBS(k).SensorParameters.sensorType == "radar")
            currTargetVisibility = 0;
        else
            currTargetVisibility = ...
                isTargetVisibleToObserver_FOV...
                (rTrg, OBS.xMat(timestep, 1:3)', ...
                dirPointing, dirSun, TRG(k).D_t, ...
                OBS.SensorParameters);
        end
        nObservedTargets = nObservedTargets + currTargetVisibility;
        if prevObsTrgCrossVisibilityVec(k)
            currTargetVisibility = ...
                currTargetVisibility * OPTIONS.valueOfPreviouslyObserved;
        end
        if currTargetVisibility > 0
            taggedTargets = [taggedTargets, k];
            targetWeights = [targetWeights, currTargetVisibility];
        end
        nWeightedTargets = nWeightedTargets + currTargetVisibility;
    end
end