function [estimationError, estimationErrorMat, traceError,...
    normalizedEstimationError] = computeEstimationError(timeVec, TRG)
    estimationError = zeros(size(timeVec));
    traceError = zeros(size(timeVec));
    normalizedEstimationError = zeros(size(timeVec));
    for j = 1:size(timeVec, 2)
        estimationError(j) = 0;
        traceError(j) = 0;
        normalizedEstimationError(j) = 0;
        for k = 1:size(TRG, 2)
            EE = TRG(k).xMatEst(j, 1:3) - TRG(k).xMat(j, 1:3);
            P = TRG(k).P(1:3,1:3,j);
            estimationErrorMat(k, j) = norm(EE);
            estimationError(j) = estimationError(j) + norm(EE);
            traceError(j) = traceError(j) + sqrt(trace(P));
            normalizedEstimationError(j) = normalizedEstimationError(j) + ...
                EE * ( P \ EE') ;
        end
        estimationError(j) = estimationError(j) / size(TRG, 2);
        traceError(j) = traceError(j) / size(TRG, 2);
        normalizedEstimationError(j) = normalizedEstimationError(j) / size(TRG, 2);
    end
end