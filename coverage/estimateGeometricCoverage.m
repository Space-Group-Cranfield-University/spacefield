function coverage = estimateGeometricCoverage(rObsMat, dirSun, nTrg, nFold, SensorParameters, D_trg, CONST, fileNames, dataFolder)
    if nargin < 9
        dataFolder = 'data';
    end
    if nargin < 8
        fileNames = "esa-master-2024-leo-"+["altitude", "declination", "size-above-10-cm"];
    end
    if nargin < 7
        CONST = initializeAstronomicalConstants();
    end
    if nargin < 6
        D_trg = 0;
    end
    if nargin < 5
        SensorParameters = getReducedSensorParameters;
    end
    if nargin < 4
        nFold = 1;
    end
    if nargin < 3
        nTrg = 1e4;
    end
    if nargin < 2
        dirSun = [-1, 0, 0]';
    end
    if ~D_trg
        integrand   = @(rTrgAugmented) ...
                    isTargetVisibleToConstellation_FOR(rTrgAugmented(1:3), ...
                    rObsMat, dirSun, rTrgAugmented(4), nFold, SensorParameters, CONST.R_E + 100);
    else
        integrand   = @(rTrgAugmented) ...
                    isTargetVisibleToConstellation_FOR(rTrgAugmented(1:3), ...
                    rObsMat, dirSun, D_trg, nFold, SensorParameters, CONST.R_E + 100);
    end
    rTrgAugmentedMat = sampleTargetPopulationFromDistribution(nTrg, fileNames, dataFolder, CONST.R_E);
    coverage = integrateMonteCarlo(rTrgAugmentedMat, integrand);
end