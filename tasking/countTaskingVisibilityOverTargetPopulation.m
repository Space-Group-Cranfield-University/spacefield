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
        currTargetVisibility = ...
            isTargetVisibleToObserver_FOV...
            (TRG(k).xMat(timestep, 1:3)', OBS.xMat(timestep, 1:3)', ...
            dirPointing, dirSun, TRG(k).D_t, ...
            OBS.SensorParameters);
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