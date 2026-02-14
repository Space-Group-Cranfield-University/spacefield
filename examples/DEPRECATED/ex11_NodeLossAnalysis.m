% This script analyses the performance loss after node loss for a given
% constellation. The visibility function neglects the Sun effects (eclipse
% and blinding) to allow generalizability of the results. On average, this
% won't affect the expected results as we are interested only in relative
% coverage loss, and not in absolute coverage estimation.

% This script replicates results shown in Fig. 22 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

%% Settings

Constants = initialiseAstronomicalConstants();
ConstellationParameters.nOrb = 15;
ConstellationParameters.nSatOrb = 27;
ConstellationParameters.h = 460;
ConstellationParameters.in = deg2rad(77.1);
ConstellationParameters.raan = 0;
ConstellationParameters.dtheta = deg2rad(8.5);
CoverageOptions.coverageCount = 3; % Triple coverage
R = Constants.R_E + 100;
deltaR = 2000;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e3; % substitute 1e4
numberLostNodes = 10; % substitute 40
lostNodesvec = 1:numberLostNodes;
nMonteCarloCoverageLoss = 5; % substitute 20
periodConstellation     = computeOrbitalPeriod(Constants.R_E + ConstellationParameters.h, ...
                        Constants.MU_E);
nTimeSteps = 4;
timeVec = linspace(0, periodConstellation, nTimeSteps + 1);
timeVec = timeVec(1:(end-1));

% Initialise Monte Carlo sample (fixed for the simulation)
rTrgMat = initialiseSampleShell(SampleOptions);

% Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);
rObsMat = convertStructToMatStatic(Constellation);

% Compute base coverage
baseCoverage = computeInstantaneousCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions);

%% Node loss analysis (about 30 s per Monte Carlo run for  timestep)

tic
coverageMat = nan(nMonteCarloCoverageLoss, numberLostNodes);
for iCoverageLoss = 1:nMonteCarloCoverageLoss
    coverageMat(iCoverageLoss, :)   = runMonteCarloCoverageLoss(numberLostNodes, ...
                                    timeVec, Constellation, rTrgMat, R, deltaR, ...
                                    CoverageOptions, Constants.MU_E);
    disp("Run: "+iCoverageLoss+" / "+ nMonteCarloCoverageLoss+", Time: "+toc)
end
meanCoverageVec = mean(coverageMat, 1);

%% Plots

figure()
hold on
for iCoverageLoss = 1:nMonteCarloCoverageLoss
    if iCoverageLoss == 1
        p_1 = stairs(lostNodesvec, coverageMat(iCoverageLoss, :) / baseCoverage * 100, "c");
    else
        stairs(lostNodesvec, coverageMat(iCoverageLoss, :) / baseCoverage * 100, "c")
    end
end
p_avg = plot(lostNodesvec, meanCoverageVec / baseCoverage * 100, "r", "LineWidth", 2);
lgd = legend([p_1,p_avg], "MC run", "mean");
lgd.Location = "northeast";
fontsize(gca, 12, 'points')
xlabel("number of lost nodes", "FontSize", 14)
ylabel("triple coverage efficiency [%]", "FontSize", 14)
xlim([1, numberLostNodes])

%% Functions

% Transform constellation structure in matrix (static constellation)
function rObsMat = convertStructToMatStatic(Constellation)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x0(1:3, 1);
    end
end

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

% Compute coverage loss with node loss
function coverageVec = runMonteCarloCoverageLoss(numberLostNodes, timeVec, Constellation, rTrgMat, R, deltaR, CoverageOptions, mu)
    coverageVec = nan(1, numberLostNodes);
    for lostNode = 1:numberLostNodes
        randomIndex = randi(size(Constellation, 2));
        Constellation(:, randomIndex) = [];
        coverageVec(lostNode) = computeCoverage(timeVec, Constellation, rTrgMat, R, deltaR, CoverageOptions, mu);
    end
end