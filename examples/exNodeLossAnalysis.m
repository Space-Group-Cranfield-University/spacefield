% This script analyses the performance loss after node loss for a given
% constellation. The visibility function neglects the Sun effects (eclipse
% and blinding) to allow generalizability of the results. On average, this
% won't affect the expected results as we are interested only in relative
% coverage loss, and not in absolute coverage estimation.

clear
clc
close all

%% Settings

Constants = initialiseAstronomicalConstants();
ConstellationParameters.nOrb = 8;
ConstellationParameters.nSatOrb = 50;
ConstellationParameters.h = 460;
ConstellationParameters.in = deg2rad(70);
ConstellationParameters.raan = 0.0625;
ConstellationParameters.dtheta = 0.0951;
CoverageOptions.coverageCount = 3; % Triple coverage
R = Constants.R_E + 100;
deltaR = 2000;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e4;
numberLostNodes = 40;
lostNodesvec = 1:numberLostNodes;
nMonteCarloCoverageLoss = 20;

% Initialise Monte Carlo sample (fixed for the simulation)
rTrgMat = initialiseSampleShell(SampleOptions);

% Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);
rObsMat = convertStructToMat(Constellation);

% Compute base coverage
baseCoverage = computeCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions);

%% Node loss analysis (about 30 s per Monte Carlo run)

tic
coverageMat = nan(nMonteCarloCoverageLoss, numberLostNodes);
for iCoverageLoss = 1:nMonteCarloCoverageLoss
    coverageMat(iCoverageLoss, :) = runMonteCarloCoverageLoss(numberLostNodes, rObsMat, rTrgMat, R, deltaR, CoverageOptions);
    disp("Run: "+iCoverageLoss+" / "+"nMonteCarloCoverageLoss, Time: "+toc)
end
meanCoverageVec = mean(coverageMat, 1);

%% Plots

figure()
hold on
for iCoverageLoss = 1:nMonteCarloCoverageLoss
    plot(lostNodesvec, coverageMat(iCoverageLoss, :) / baseCoverage * 100, "c")
end
plot(lostNodesvec, meanCoverageVec / baseCoverage * 100, "r")
xlabel("Number of lost nodes")
ylabel("Coverage loss [%]")

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
function flag = isVisibleNeglectingSun(rObs, rTrg, R, deltaR)
    flag = ~isWithinRange(rObs, rTrg, deltaR);
    if flag == 0
        flag = isObstructed(rObs, rTrg, R);
    end
    flag = ~flag;
end

% Define instantaneous coverage function
function coverage = computeCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions)
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

% Compute coverage loss with node loss
function coverageVec = runMonteCarloCoverageLoss(numberLostNodes, rObsMat, rTrgMat, R, deltaR, CoverageOptions)
    coverageVec = nan(1, numberLostNodes);
    for lostNode = 1:numberLostNodes
        randomIndex = randi(size(rObsMat, 1));
        rObsMat(randomIndex, :) = [];
        coverageVec(lostNode) = computeCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions);
    end
end