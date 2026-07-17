function [TRG, nWellTracked, nWellTrackedVec, isWellTrackedMat] = ...
            computeNumberOfWellTrackedTargets...
            (TRG, targetVisibilityMat, bestAvailableSensor, deltaT, minTrackLength)
    if nargin < 5
        minTrackLength = 4 * deltaT;
    end
    TRG = extractTracksFromVisibilityMatrix(TRG, targetVisibilityMat);
    isWellTrackedMat = zeros(size(targetVisibilityMat));
    for k = 1:size(TRG, 2)
        T = computeOrbitalPeriod(TRG(k).a);
        timeAheadOptical = floor(0.5 * T / deltaT);
        timeAheadRadar = floor(0.5 * T / deltaT);
        for i = 1:size(TRG(k).track, 1)
            firstIndex = TRG(k).track(i, 1);
            lastIndex = TRG(k).track(i, 2);
            trackLength = ( 1 + lastIndex - firstIndex ) * deltaT;
            lastIndexOptical = lastIndex;
            lastIndexRadar = -1e9;
            for j = firstIndex:lastIndex
                bestSensor = bestAvailableSensor(k, firstIndex);
                if bestSensor == "radar"
                    lastIndexRadar = j;
                end
            end
            if trackLength >= minTrackLength
                lastIndexSupport = max(lastIndexOptical + timeAheadOptical, lastIndexRadar + timeAheadRadar);
                lastIndexSupport = min(lastIndexSupport, size(targetVisibilityMat, 2));
                isWellTrackedMat(k, firstIndex:lastIndexSupport) = 1;
            end
        end
    end
    T = computeOrbitalPeriod(TRG(k).a);
    timestep = floor(T / deltaT);
    nWellTrackedVec = sum(isWellTrackedMat, 1);
    nWellTracked = mean(nWellTrackedVec(timestep:end));
end