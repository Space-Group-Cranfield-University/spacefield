function bestGndLocationMat = searchOptimalRadarNetwork(N, OPTIONS)
    if nargin < 2
        OPTIONS = getStandardRadarSearchOptions;
    end
    nRuns = OPTIONS.nRuns;
    RadarParameters = OPTIONS.RadarParameters;
    dirSun = OPTIONS.dirSun;
    nFold = OPTIONS.nFold;
    nTrg = OPTIONS.nTrg;
    
    bestCoverage = 0;
    bestGndLocationMat = 0;
    for k = 1:nRuns
        gndLocationMat = getRandomGroundLocations(N);
        [~, rObsMat] = initializeGroundNetwork(gndLocationMat, RadarParameters);
        coverage = estimateGeometricCoverage(rObsMat, dirSun, nTrg, nFold, RadarParameters);
        if coverage > bestCoverage
            bestCoverage = coverage;
            bestGndLocationMat = gndLocationMat;
        end
    end
end