function [ConstellationParameters, bestCoverageValue] = ...
    searchOptimalConstellation(N, OPTIONS)
    if nargin < 2
        OPTIONS = getStandardConstellationSearchOptions;
    end
    % Parameter extraction
    dN = OPTIONS.dN;
    nRuns = OPTIONS.nRuns;
    nTests = OPTIONS.nTests;
    nFold = OPTIONS.nFold;
    dirSun0 = OPTIONS.dirSun0;
    n_t_season = OPTIONS.n_t_season;
    n_t_day = OPTIONS.n_t_day;
    isSeasonal = OPTIONS.isSeasonal;
    isDaily = OPTIONS.isDaily;
    nTrg = OPTIONS.nTrg;
    SensorParameters = OPTIONS.SensorParameters;
    
    % Search
    bestCoverageMat = zeros(nTests, nRuns);
    meanCoverageVec = zeros(1, nRuns);
    for k = 1:nRuns
        for j = 1:nTests
            Parameters = sampleRandomConstellationParameters([N-dN,N]);
            OBS = initializeWalkerConstellation(Parameters);
            if isSeasonal
                coverage = getSeasonalCoverage(OBS, dirSun0, n_t_season, n_t_day, nFold);
            elseif isDaily
                timeVecDay = linspace(0, OBS(1).T, n_t_day);
                coverage = getDailyCoverage(timeVecDay, OBS, dirSun0, nFold, nTrg);
            else
                rObsMat = getConstellationPositionMatrix(OBS);
                coverage = estimateGeometricCoverage(rObsMat, dirSun0, nTrg, nFold, SensorParameters);
            end
            if k == 1 || coverage > bestCoverageMat(j, k-1)
                bestCoverageMat(j, k) = coverage;
                bestParameters(j) = Parameters;
            else
                bestCoverageMat(j, k) = bestCoverageMat(j, k-1);
            end
        end
    %disp("Run: "+string(k)+" / "+string(nRuns))
    meanCoverageVec(k) = mean(bestCoverageMat(:, k));
    end
    indexBestCoverage = find(bestCoverageMat(:, end) == max(bestCoverageMat(:, end)));
    bestCoverageValue = bestCoverageMat(indexBestCoverage, end);
    ConstellationParameters = bestParameters(indexBestCoverage);
end