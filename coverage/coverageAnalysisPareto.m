function RESULTS = coverageAnalysisPareto(nConstellations, nFold, nSatBounds, inBounds, hBounds, excludeSingleOrbits, n_t_season, n_t_day)
    if nargin < 8
        n_t_day = 3;
    end
    if nargin < 7
        n_t_season = 4;
    end
    if nargin < 6
        excludeSingleOrbits = 1;
    end
    if nargin < 5
        hBounds = 460;
    end
    if nargin < 4
        inBounds = [pi/2 - pi/6, pi/2];
    end
    if nargin < 3
        nSatBounds = [2 100];
    end
    if nargin < 2
        nFold = 1;
    end
    dirSun0 = [-1 0 0]';
    nSat = zeros(1, nConstellations);
    in = zeros(1, nConstellations);
    h = zeros(1, nConstellations);
    coverage = zeros(1, nConstellations);
    stdCoverage = zeros(1, nConstellations);
    tic
    for k = 1:nConstellations
        Parameters = sampleRandomConstellationParameters(nSatBounds, inBounds, hBounds, excludeSingleOrbits);
        OBS = initializeWalkerConstellation(Parameters);
        nSat(k) = Parameters.nOrb * Parameters.nSatOrb;
        in(k) = rad2deg(Parameters.in);
        h(k) = Parameters.h;
        [coverage(k), varCoverage] = getSeasonalCoverage(OBS, dirSun0, n_t_season, n_t_day, nFold);
        stdCoverage(k) = sqrt(varCoverage);
        currentTime = toc;
        if ~mod(k, 20)
            disp(string(k)+" / "+string(nConstellations)+", time elapsed [s]: "+string(currentTime))
        end
    end
    [nSatPareto, coveragePareto] = estimateParetoFrontier(nSat, coverage);
    RESULTS.nSat = nSat;
    RESULTS.in = in;
    RESULTS.h = h;
    RESULTS.coverage = coverage;
    RESULTS.stdCoverage = stdCoverage;
    RESULTS.coveragePareto = coveragePareto;
    RESULTS.nSatPareto = nSatPareto;
end