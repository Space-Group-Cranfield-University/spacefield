function [nTrackedTargets, nTrackedTargetsVec] = computeNumberOfTargetsTrackedOnce(targetVisibilityMat, nMeasurements)
    if nargin < 2
        nMeasurements = 1;
    end
    for j = 1:size(targetVisibilityMat, 2)
        nTrackedTargetsVec(j) = countTrackedTargets(targetVisibilityMat(:, 1:j), nMeasurements);
    end
    nTrackedTargets = mean(nTrackedTargetsVec);
end

function nTrackedTargets = countTrackedTargets(targetVisibilityMat, nMeasurements)
    targetVisibilityMat = (targetVisibilityMat >= nMeasurements);
    targetVisibilityVec = sum(targetVisibilityMat, 2);
    targetVisibilityVec = (targetVisibilityVec >= 1);
    nTrackedTargets = sum(targetVisibilityVec);
end