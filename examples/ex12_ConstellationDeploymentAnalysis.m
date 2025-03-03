% This script compares different deployment strategies for a given
% constellation. The most effective deployment is the one that achieves a
% working coverage faster. The visibility function neglects the Sun effects 
% (eclipse and blinding) to allow generalizability of the results. On 
% average, this won't affect the expected results as we are interested only
% in the relative coverage loss, and not in absolute coverage estimation.
% We assume a single launch brings 10 satellites to orbit (40 launches
% necessary to complete the constellation). Since the series of 10 deployed
% satellites is not evenly distributed on the orbit, a time averaging is
% needed for computing the actual average coverage of the constellation.

% This script replicates results shown in Fig. 24 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

%% Settings

Constants = initialiseAstronomicalConstants();
ConstellationParameters.nOrb = 15;
ConstellationParameters.nSatOrb = 27;
ConstellationParameters.h = 460;
ConstellationParameters.in = deg2rad(77.1);
ConstellationParameters.raan = deg2rad(0);
ConstellationParameters.dtheta = deg2rad(8.5);
CoverageOptions.coverageCount = 3; % Triple coverage
R = Constants.R_E + 100;
deltaR = 2000;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e3; % substitute 1e4
nSatLaunch = 9;
nSat = ConstellationParameters.nOrb * ConstellationParameters.nSatOrb;
nLaunches = nSat / nSatLaunch;
nSeriesOrb = nLaunches / ConstellationParameters.nOrb;
launchesVec = 1:nLaunches;
periodConstellation     = computeOrbitalPeriod(Constants.R_E + ConstellationParameters.h, ...
                        Constants.MU_E);
nTimeSteps = 4;
timeVec = linspace(0, periodConstellation, nTimeSteps + 1);
timeVec = timeVec(1:(end-1));
strategyVec = ["fullOrbitFirst", "subsequentOrbitsWithSkip",...
            "randomSeries"];
for k = 1:size(strategyVec, 2)
    legendVec(k) = "strategy "+string(k);
end

legendVec = ["fullOrbitFirst", "consecutiveOrbits", "random"];
lastLaunch = ConstellationParameters.nOrb * ConstellationParameters.nSatOrb / nSatLaunch;

%% Deployment strategy analysis

% Initialise Monte Carlo sample
rTrgMat = initialiseSampleShell(SampleOptions);

% Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);

coverageMat = nan(size(strategyVec, 2), nLaunches);
tic
for k = 1:size(strategyVec, 2)
    strategy = strategyVec(k);
    coverageMat(k, :)   = deployConstellation(nLaunches, timeVec, Constellation, ...
                        rTrgMat, R, deltaR, CoverageOptions, Constants.MU_E, ...
                        ConstellationParameters, nSatLaunch, nSeriesOrb, strategy);
    disp("Completed: "+k+" / "+size(strategyVec, 2)+", time: "+toc)
end

%% Plots

figure()
for k = 1:size(strategyVec, 2)
    plot(launchesVec, coverageMat(k, :) * 100 / coverageMat(k, end), "LineWidth", 2)
    hold on
end
fontsize(gca, 12, 'points')
xlabel("number of launches", "FontSize", 14)
ylabel("triple coverage efficiency [%]", "FontSize", 14) %  cov_{3,deployed)}/ cov_{3,full}
title("Constellation deployment analysis")
lgd = legend(legendVec);
lgd.Location = "southeast";
xlim([1, lastLaunch])

%% Functions

% Transform constellation structure in matrix
function rObsMat = convertStructToMat(Constellation, timeStep)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x(timeStep, 1:3)';
    end
end

% Define visibility function
function flag = isVisibleNeglectingSun(rObs, rTrg, R, deltaR)
    flag = ~isWithinRange(rObs, rTrg, deltaR);
    if flag == 0
        flag = ~isAboveTheHorizon(rObs, rTrg, R);
    end
    flag = ~flag;
end

% Define instantaneous coverage function
function coverage = computeInstantaneousCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions)
    nObs = size(rObsMat, 1);
    nTrg = size(rTrgMat, 1);
    coverage = 0;
    for iTrg = 1:nTrg
        rTrg = rTrgMat(iTrg,:)';
        count = 0;
        iObs = 1;
        while count < CoverageOptions.coverageCount && iObs < nObs
            rObs = rObsMat(iObs,:)';
            count = count + isVisibleNeglectingSun(rObs, rTrg, R, deltaR);
            iObs = iObs + 1;
        end
        coverage = coverage + (count == CoverageOptions.coverageCount);
    end
    coverage = coverage / nTrg;
end

% Define average coverage function
function coverage = computeCoverage(timeVec, Constellation, rTrgMat, R, deltaR, CoverageOptions, mu)
    Constellation = propagateConstellation(timeVec, Constellation, mu);
    totalCoverage = 0;
    for timeStep = 1:size(timeVec, 2)
        rObsMat = convertStructToMat(Constellation, timeStep);
        totalCoverage = totalCoverage + computeInstantaneousCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions);
    end
    coverage = totalCoverage / size(timeVec, 2);
end

% Function for getting deployment coverage vector
function coverageVec    = deployConstellation(nLaunches, timeVec, Constellation, ...
                        rTrgMat, R, deltaR, CoverageOptions, mu, ...
                        ConstellationParameters, nSatLaunch, nSeriesOrb, strategy)
    CurrentConstellation = [];
    coverageVec = nan(1, nLaunches);
    if strcmp(strategy, "randomSeries")
        launchVec = 1:nLaunches;
    else
        launchVec = 1:size(Constellation, 2);
    end
    for iLaunches = 1:nLaunches
        if strcmp(strategy, "fullOrbitFirst")
            seriesLaunched = 1 + (iLaunches - 1)*nSatLaunch : iLaunches*nSatLaunch;
        end
        if strcmp(strategy, "subsequentOrbits")
            currentOrbit = mod(iLaunches, ConstellationParameters.nOrb);
            if currentOrbit == 0
                currentOrbit = ConstellationParameters.nOrb;
            end
            currentSeries = ceil(iLaunches / ConstellationParameters.nOrb);
            seriesLaunched  = (currentOrbit - 1)*ConstellationParameters.nSatOrb + ...
                            (currentSeries - 1)*nSatLaunch + 1 : ...
                            (currentOrbit - 1)*ConstellationParameters.nSatOrb + ...
                            currentSeries*nSatLaunch;
        end
        if strcmp(strategy, "subsequentOrbitsWithSkip")
            currentOrbit = mod(iLaunches, ConstellationParameters.nOrb);
            if currentOrbit == 0
                currentOrbit = ConstellationParameters.nOrb;
            end
            currentSeries = ceil(iLaunches / ConstellationParameters.nOrb);
            skipOrbit = mod(mod(currentOrbit, nSeriesOrb) - 1, nSeriesOrb); % Correct    
            seriesLaunched  = (currentOrbit - 1)*ConstellationParameters.nSatOrb + ...
                            mod((currentSeries + skipOrbit - 1), nSeriesOrb)*nSatLaunch + 1 : ...
                            (currentOrbit - 1)*ConstellationParameters.nSatOrb + ...
                            (mod((currentSeries + skipOrbit - 1), nSeriesOrb) + 1)*nSatLaunch;
        end
        if strcmp(strategy, "randomSeries")
            currentLaunch = launchVec(randi(size(launchVec, 2)));
            launchVec(launchVec == currentLaunch) = [];
            seriesLaunched = (currentLaunch - 1) * nSatLaunch + 1 : currentLaunch * nSatLaunch;
        end
        if strcmp(strategy, "randomSats")
            seriesLaunched = zeros(1, nSatLaunch);
            for iSeries = 1:nSatLaunch
                seriesLaunched(iSeries) = launchVec(randi(size(launchVec, 2)));
                launchVec(launchVec == seriesLaunched(iSeries)) = [];
            end
        end
        CurrentConstellation = [CurrentConstellation, Constellation(seriesLaunched)];
        coverageVec(iLaunches)  = computeCoverage(timeVec, CurrentConstellation, ...
                                rTrgMat, R, deltaR, CoverageOptions, mu);
    end
end