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
            trackLength = ( lastIndex - firstIndex ) * deltaT;
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


function [TRG, nWellTracked, nWellTrackedVec, isWellTrackedMat] = ...
            computeNumberOfWellTrackedTargetsBackUp...
            (TRG, targetVisibilityMat, bestAvailableSensor, deltaT, maxOrbitalFraction, minTrackLength)
    if nargin < 5
        minTrackLength = 4 * deltaT;
    end
    if nargin < 4
        maxOrbitalFraction = 1/2;
    end
    TRG = extractTracksFromVisibilityMatrix(TRG, targetVisibilityMat);
    isWellTrackedMat = zeros(size(targetVisibilityMat));
    for k = 1:size(TRG, 2)
        T = computeOrbitalPeriod(TRG(k).a);
        for j = 1:size(targetVisibilityMat, 2)
            bestSensor = bestAvailableSensor(k, j);
            lastMeasurementTimestep = -1e9;
            for i = 1:size(TRG(k).track, 1)
                trackLength = ( TRG(k).track(i, 2) - TRG(k).track(i, 1) ) * deltaT;
                if trackLength >= minTrackLength
                    if j >= TRG(k).track(i, 1) && j <= TRG(k).track(i, 2)
                        isWellTrackedMat(k, j) = 1;
                    elseif j > TRG(k).track(i, 2)
                        lastMeasurementTimestep = TRG(k).track(i, 2);
                    end
                end
            end
            timeSinceLastMeasurement = ...
                (j - lastMeasurementTimestep) * deltaT;
            if timeSinceLastMeasurement < ( maxOrbitalFraction * T )
                isWellTrackedMat(k, j) = 1;
            end
        end
    end
    nWellTrackedVec = sum(isWellTrackedMat, 1);
    nWellTracked = mean(nWellTrackedVec);
end