function [meanWellTrackedEstimationError, meanWellTrackedEstimationErrorVec] = ...
    computeEstimationErrorOverWellTrackedTargets(estimationErrorMat, isWellTrackedMat)
    meanWellTrackedEstimationErrorVec = nan(1, size(isWellTrackedMat, 2));
    meanWellTrackedEstimationError = 0;
    for j = 1:size(isWellTrackedMat, 2)
        %meanWellTrackedEstimationErrorVec(j) = 0;
        nWellTrackedTargets = 0;
        for k = 1:size(isWellTrackedMat, 1)
            if isWellTrackedMat(k, j)
                if isnan(meanWellTrackedEstimationErrorVec(j))
                    meanWellTrackedEstimationErrorVec(j) = ...
                        estimationErrorMat(k, j);
                    nWellTrackedTargets = 1;
                else
                    nWellTrackedTargets = nWellTrackedTargets + 1;
                    meanWellTrackedEstimationErrorVec(j) = ...
                        ( (nWellTrackedTargets - 1) * ...
                        meanWellTrackedEstimationErrorVec(j) + ...
                        estimationErrorMat(k, j) ) / nWellTrackedTargets;
                end
            end
        end
    end
end