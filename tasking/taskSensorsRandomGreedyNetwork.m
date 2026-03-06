function OBS = taskSensorsRandomGreedyNetwork...
        (timestep, OBS, TRG, prevObsTrgCrossVisibilityMat, dirSun, OPTIONS)
    
    bestObservedTargets = 0;
    bestPointingDirSet = zeros(3, size(OBS, 2));
    bestAzElSet = zeros(2, size(OBS, 2));
    for j = 1:OPTIONS.nSamples
        observedTargetVec = zeros(1, size(TRG, 2));
        currentPointingDirectionSet = [];
        currentAzElSet = [];
        for k = 1:size(OBS, 2)
            prevPointingDir = OBS(k).pointingMat(timestep - 1, :)';
            prevAzEl = convertPositionToRaDec(prevPointingDir);
            if abs(prevAzEl - 1) > 1e-6
                warning("unit vector not unitary")
            end
            prevAzEl = prevAzEl(2:3);
            if j == 1
                deltaAzEl = [0;0];
            elseif OPTIONS.constrained
                deltaAzEl = sampleAzElGivenMaxSlewAngle(prevAzEl, OPTIONS, 1);
            else
                deltaAzEl = ( 2*rand(2, 1) - 1 ) * OPTIONS.maxSlewAngle;
            end
            newAzEl = prevAzEl + deltaAzEl;
            newPointingDir = convertRaDecToPosition([1; newAzEl]);
            currentAzElSet = [currentAzElSet, newAzEl];
            currentPointingDirectionSet = [currentPointingDirectionSet, newPointingDir];
            OBS(k).pointingMat(timestep, :) = newPointingDir';
            [~, ~, taggedTargets, targetWeights] = ...
                countTaskingVisibilityOverTargetPopulation...
                (timestep, OBS(k), TRG, prevObsTrgCrossVisibilityMat(k, :), ...
                dirSun, OPTIONS);
            observedTargetTemp = observedTargetVec(taggedTargets);
            observedTargetTemp(observedTargetTemp < 2) = targetWeights;
            observedTargetVec(taggedTargets) = observedTargetTemp;
        end
        nObservedTargets = sum(observedTargetVec > 0);
        % disp("Sample: "+string(j)+", Current observed: "+string(nObservedTargets)+", Best observed: "+string(bestObservedTargets))
        if j == 1 || nObservedTargets > bestObservedTargets
            bestObservedTargets = nObservedTargets;
            bestPointingDirSet = currentPointingDirectionSet;
            bestAzElSet = currentAzElSet;
        end
    end
    for k = 1:size(OBS, 2)
        OBS(k).pointingMat(timestep, :) = bestPointingDirSet(:, k)';
        OBS(k).pointingAzElMat(timestep, :) = bestAzElSet(:, k)';
    end
end