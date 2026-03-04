function [completeness, completenessVec] = computeCompleteness(targetVisibilityMat, completenessFold)
    if nargin < 2
        completenessFold = 1;
    end
    completeness = computeCompletenessAtFixedTimestep...
        (targetVisibilityMat, completenessFold);
    for j = 1:size(targetVisibilityMat, 2)
        completenessVec(j) = computeCompletenessAtFixedTimestep...
            (targetVisibilityMat, completenessFold, j);
    end
end

function completeness = ...
    computeCompletenessAtFixedTimestep(targetVisibilityMat, completenessFold, timestep)
    if nargin < 2
        completenessFold = 1;
    end
    if nargin < 3
        timestep = size(targetVisibilityMat, 2);
    end
    completeness = ...
        sum(sum(targetVisibilityMat(:, 1:timestep) >= completenessFold)) / ...
        numel(targetVisibilityMat(:, 1:timestep));
end