function OBS = taskSensorsRandomGreedySingleSensor...
    (timestep, OBS, TRG, prevObsTrgCrossVisibilityMat, dirSun, OPTIONS)

    for k = 1:size(OBS, 2)
        prevPointingDir = OBS(k).pointingMat(timestep - 1, :)';
        prevAzEl = convertPositionToRaDec(prevPointingDir);
        if abs(prevAzEl - 1) > 1e-6
            warning("unit vector not unitary")
        end
        prevAzEl = prevAzEl(2:3);
        bestWeightedTargets = 0;
        bestPointingDir = prevPointingDir;
        bestAzEl = prevAzEl;
        if OPTIONS.constrained
            deltaAzElVec = sampleAzElGivenMaxSlewAngle(prevAzEl, OPTIONS);
        else
            %deltaAzElVec = rand([2, OPTIONS.nSamples]) .* deg2rad(90) - [0; deg2rad(45)];
            deltaAzElVec = (2*rand(2,OPTIONS.nSamples)-1) * OPTIONS.maxSlewAngle;
        end
        for j = 1:OPTIONS.nSamples
            % We first check static pointing to avoid slewing when 
            % unnecessary.
            if j > 1
                deltaAzEl = deltaAzElVec(:,j);
            else
                deltaAzEl = [0;0];
            end
            newAzEl = prevAzEl + deltaAzEl;
            newPointingDir = convertRaDecToPosition([1; newAzEl]);
            OBS(k).pointingMat(timestep, :) = newPointingDir';
            [~, nWeightedTargets, ~] = ...
                countTaskingVisibilityOverTargetPopulation...
                (timestep, OBS(k), TRG, prevObsTrgCrossVisibilityMat(k, :), ...
                dirSun, OPTIONS);
            % We select the pointing direction that maximizes the number of
            % observed targets for the single sensor, weighted over target
            % priority (targets observed at a previous time step have a
            % higher priority).
            if nWeightedTargets > bestWeightedTargets
                bestWeightedTargets = nWeightedTargets;
                bestPointingDir = newPointingDir;
                bestAzEl = newAzEl;
            end
        end
        OBS(k).pointingMat(timestep, :) = bestPointingDir';
        OBS(k).pointingAzElMat(timestep, :) = bestAzEl';
    end
end