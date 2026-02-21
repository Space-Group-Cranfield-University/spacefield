function [meanCoverage, varCoverage, dailyCoverageVec] = getDailyCoverage(timeVecDay, OBS, dirSun, nFold)
    OBS = propagateConstellation(timeVecDay, OBS);
    dailyCoverageVec = zeros(1, size(timeVecDay, 2));
    for k = 1:size(timeVecDay, 2)
        rObsMat = getConstellationPositionMatrix(OBS, k);
        dailyCoverageVec(k) = estimateGeometricCoverage(rObsMat, dirSun, 1e4, nFold);
        if ~mod(k, 100)
            disp(string(k)+" / "+string(size(timeVecDay, 2)));
        end
    end
    meanCoverage = mean(dailyCoverageVec);
    varCoverage = var(dailyCoverageVec);
end