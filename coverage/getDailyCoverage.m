function [meanCoverage, varCoverage, dailyCoverageVec] = getDailyCoverage(timeVecDay, OBS, dirSun, nFold, nSamples)
    if nargin < 4
        nFold = 1;
    end
    if nargin < 3
        dirSun = [-1, 0, 0]';
    end
    if nargin < 5
        nSamples = 1e4;
    end
    if ~isfield(OBS, "xMat") && OBS(1).SensorParameters.sensorType == "optical"
        OBS = propagateConstellation(timeVecDay, OBS);
    elseif ~isfield(OBS, "xMat") && OBS(1).SensorParameters.sensorType == "radar"
        OBS = propagateGroundNetwork(timeVecDay, OBS);
    end
    dailyCoverageVec = zeros(1, size(timeVecDay, 2));
    for k = 1:size(timeVecDay, 2)
        rObsMat = getConstellationPositionMatrix(OBS, k);
        dailyCoverageVec(k) = estimateGeometricCoverage(rObsMat, dirSun, nSamples, nFold, OBS(1).SensorParameters);
        if ~mod(k, 100)
            disp(string(k)+" / "+string(size(timeVecDay, 2)));
        end
    end
    meanCoverage = mean(dailyCoverageVec);
    varCoverage = var(dailyCoverageVec);
end