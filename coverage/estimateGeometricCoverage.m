function coverage = estimateGeometricCoverage(nTrg, rObsMat, dirSun, nFold, V_lim, CONST, fileNames, dataFolder)
    if nargin < 8
        dataFolder = 'data';
    end
    if nargin < 7
        fileNames = "esa-master-2024-leo-"+["altitude", "declination", "size-above-10-cm"];
    end
    if nargin < 6
        CONST = initializeAstronomicalConstants();
    end
    if nargin < 5
        V_lim = 16;
    end
    if nargin < 4
        nFold = 1;
    end
    integrand   = @(rTrgAugmented) ...
                isTargetVisibleToConstellation_FOR(rTrgAugmented(1:3), ...
                rObsMat, dirSun, rTrgAugmented(4), nFold, V_lim, CONST.R_E + 100);
    rTrgAugmentedMat = sampleTargetPopulationFromDistribution(nTrg, fileNames, dataFolder, CONST.R_E);
    coverage = integrateMonteCarlo(rTrgAugmentedMat, integrand);
end