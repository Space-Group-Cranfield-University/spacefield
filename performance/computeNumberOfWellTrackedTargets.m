function [TRG, nWellTracked, nWellTrackedVec, isWellTrackedMat] = ...
            computeNumberOfWellTrackedTargets(TRG, targetVisibilityMat, deltaT, orbitalFraction)
    if nargin < 4
        orbitalFraction = 1/2;
    end
    TRG = extractTracksFromVisibilityMatrix(TRG, targetVisibilityMat);
    isWellTrackedMat = zeros(size(targetVisibilityMat));
    for k = 1:size(TRG, 2)
        T = computeOrbitalPeriod(TRG(k).a);
        for j = 1:size(targetVisibilityMat, 2)
            lastMeasurementTimestep = -1e9;
            for i = 1:size(TRG(k).track, 1)
                if j >= TRG(k).track(i, 1) && j <= TRG(k).track(i, 2)
                    isWellTrackedMat(k, j) = 1;
                elseif j > TRG(k).track(i, 2)
                    lastMeasurementTimestep = TRG(k).track(i, 2);
                end
            end
            timeSinceLastMeasurement = ...
                (j - lastMeasurementTimestep) * deltaT;
            if timeSinceLastMeasurement < ( orbitalFraction * T )
                isWellTrackedMat(k, j) = 1;
            end
        end
    end
    nWellTrackedVec = sum(isWellTrackedMat, 1);
    nWellTracked = mean(nWellTrackedVec);
end