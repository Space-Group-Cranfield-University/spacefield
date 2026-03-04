function [meanEstimationErrorUnderTracking, meanEstimationErrorUnderTrackingVec, ...
        medianEstimationErrorUnderTracking, medianEstimationErrorUnderTrackingVec] = ...
        computeEstimationErrorUnderTracking...
        (estimationErrorMat, targetVisibilityMat, nFold)
    if nargin < 3
        nFold = 1;
    end
    targetVisibilityMat = (targetVisibilityMat >= nFold);
    estimationErrorMat = estimationErrorMat .* targetVisibilityMat;
    EE_list = [];
    for j =1:size(estimationErrorMat, 2)
        EE_vec = [];
        for k = 1:size(estimationErrorMat, 1)
            EE = estimationErrorMat(k, j);
            if EE > 0
                EE_vec = [EE_vec, EE];
            end
        end
        meanEstimationErrorUnderTrackingVec(j) = mean(EE_vec);
        medianEstimationErrorUnderTrackingVec(j) = median(EE_vec);
        if ~isempty(EE_vec)
            EE_list = [EE_list, EE_vec];
        end
    end
    meanEstimationErrorUnderTracking = mean(EE_list);
    medianEstimationErrorUnderTracking = median(EE_list);
end

function [meanEstimationErrorUnderTracking, ...
    medianEstimationErrorUnderTracking] = ...
    computeEstimationErrorScalar(estimationErrorMat)
    
end