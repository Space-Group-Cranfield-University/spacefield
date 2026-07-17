function [minConstellationDistanceMat, taggedSatMat] = ...
    computeMinimumDistanceBetweenConstellationSatellites(OBS)
    minConstellationDistanceMat = 1e5*ones(size(OBS(1).xMat, 1), size(OBS, 2));
    taggedSatMat = 1e5*ones(size(OBS(1).xMat, 1), size(OBS, 2));
    for k = 1:size(OBS, 2)
        rPrimaryMat = OBS(k).xMat(:, 1:3);
        minConstellationDistanceVec = 1e5*ones(size(OBS(k).xMat, 1), 1);
        for j = 1:(size(OBS, 2)-1)
            index = j;
            if j >= k
                index = index + 1;
            end
            rSecondaryMat = OBS(index).xMat(:, 1:3);
            deltaMat = rPrimaryMat - rSecondaryMat;
            for l = 1:size(deltaMat, 1)
                deltaRcurrent = norm(deltaMat(l, :));
                if deltaRcurrent < minConstellationDistanceVec(l)
                    minConstellationDistanceVec(l) = deltaRcurrent;
                    taggedSatMat(l, k) = index;
                end
            end
        end
        minConstellationDistanceMat(:, k) = minConstellationDistanceVec;
    end
end