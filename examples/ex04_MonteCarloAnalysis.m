% This script analyses the performance of Monte Carlo integration for the
% coverage computation of a given optical surveillance constellation.

% This script replicates results shown in Fig. 12 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

%% Settings

nMCRuns = 10; % Number of Monte Carlo runs per nMC. Substitute 50
nMCVec = ceil(10.^(3:0.5:4)); % substitute ceil(10.^(3:0.25:5))
Constants = initialiseAstronomicalConstants;
ConstellationParameters.nOrb = 8;
ConstellationParameters.nSatOrb = 50;
ConstellationParameters.h = 460;
ConstellationParameters.in = deg2rad(70);
ConstellationParameters.raan = 0;
ConstellationParameters.dtheta = deg2rad(5);
CoverageOptions.coverageCount = 3; % Triple coverage
R = Constants.R_E + 100;
deltaR = 2000;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
dirSun = [-1 0 0]';

%% Validate Monte Carlo method

% Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);
rObsMat = convertStructToMat(Constellation);

% Monte Carlo integration
tic
meanCoverageVec = nan(1, size(nMCVec, 2));
varCoverageVec = nan(1, size(nMCVec, 2));
for kMC = 1:size(nMCVec, 2)
    SampleOptions.N = nMCVec(kMC);
    coverageVec = nan(1, nMCRuns);
    for kRun = 1:nMCRuns
        % Initialise Monte Carlo sample
        rTrgMat = initialiseSampleShell(SampleOptions);

        % Compute coverage
        coverageVec(kRun) = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
    end
    meanCoverageVec(kMC) = mean(coverageVec);
    varCoverageVec(kMC) = var(coverageVec);
    disp("nMC: "+nMCVec(kMC)+", time: "+toc)
end

%% Plot

figure()
loglog(nMCVec, varCoverageVec, "LineWidth", 2)
hold on
loglog(nMCVec, 1./nMCVec, "LineWidth", 2)
legend(["MC Sample Variance","1/n_{MC}"])
fontsize(gca, 12, 'points')
xlabel("number of MC samples", "FontSize", 14)
ylabel("triple coverage variance [-]", "FontSize", 14)
title("Triple coverage variance against Number of MC samples", "FontSize", 12)

figure()
loglog(nMCVec, 100*sqrt(varCoverageVec), "LineWidth", 2)
hold on
loglog(nMCVec, 100*sqrt(1./nMCVec), "LineWidth", 2)
legend(["MC Sample Std. Dev.","n_{MC}^{-1/2}"])
fontsize(gca, 12, 'points')
xlabel("number of MC samples", "FontSize", 14)
ylabel("triple coverage sample std. dev. [%]", "FontSize", 14)
title("Triple coverage standard deviation against Number of MC samples", "FontSize", 12)

figure()
loglog(nMCVec, 100*sqrt(meanCoverageVec), "LineWidth", 2)
fontsize(gca, 12, 'points')
xlabel("number of MC samples", "FontSize", 14)
ylabel("triple coverage sample mean [%]", "FontSize", 14)
title("Average triple coverage against Number of MC samples", "FontSize", 12)

figure()

subplot(2,1,1)
loglog(nMCVec, 100*sqrt(meanCoverageVec), "LineWidth", 2)
fontsize(gca, 12, 'points')
xlabel("number of MC samples", "FontSize", 14)
ylabel("triple coverage sample mean [%]", "FontSize", 14)
title("Average triple coverage against Number of MC samples", "FontSize", 12)

subplot(2,1,2)
loglog(nMCVec, 100*sqrt(varCoverageVec), "LineWidth", 2)
hold on
loglog(nMCVec, 100*sqrt(1./nMCVec), "LineWidth", 2)
legend(["MC Sample Std. Dev.","n_{MC}^{-1/2}"])
fontsize(gca, 12, 'points')
xlabel("number of MC samples", "FontSize", 14)
ylabel("triple coverage sample std. dev. [%]", "FontSize", 14)
title("Triple coverage standard deviation against Number of MC samples", "FontSize", 12)

%% Functions

% Transform constellation structure in matrix
function rObsMat = convertStructToMat(Constellation)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x0(1:3,1);
    end
end

% Define visibility function
function flag = isVisibleToObserver(rObs, rTrg, dirSun, R, deltaR)
    flag    = ~( ( ~isWithinRange(rObs, rTrg, deltaR) ) || ... 
            isScatteringWeak(rObs, rTrg, dirSun) || isEclipsed(rTrg, dirSun, R) || ...
            ~isAboveTheHorizon(rObs, rTrg, R) );
end

% Define aggregate visibility function
function flag = isVisibleToConstellation(rTrg, rObsMat, nObs, dirSun, R, deltaR, CoverageOptions)
    iObs = 1;
    count = 0;
    while count < CoverageOptions.coverageCount && iObs < nObs
        rObs = rObsMat(iObs,:)';
        count = count + isVisibleToObserver(rObs, rTrg, dirSun, R, deltaR);
        iObs = iObs + 1;
    end
    flag = (count == CoverageOptions.coverageCount);
end

% Define instantaneous coverage function(s)
function coverage = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions)
    nObs = size(rObsMat, 1);
    integrandFunction = @(r) isVisibleToConstellation(r, rObsMat, nObs, dirSun, R, deltaR, CoverageOptions);
    coverage = integrateMonteCarlo(rTrgMat, integrandFunction);
end

function coverage = computeCoverage_2(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions)
    nObs = size(rObsMat, 1);
    nTrg = size(rTrgMat, 1);
    coverage = 0;
    for iTrg = 1:nTrg
        rTrg = rTrgMat(iTrg,:)';
        count = 0;
        iObs = 1;
        while count < CoverageOptions.coverageCount && iObs < nObs
            rObs = rObsMat(iObs,:)';
            count = count + isVisibleToObserver(rObs, rTrg, dirSun, R, deltaR);
            iObs = iObs + 1;
        end
        coverage = coverage + (count == CoverageOptions.coverageCount);
    end
    coverage = coverage / nTrg;
end