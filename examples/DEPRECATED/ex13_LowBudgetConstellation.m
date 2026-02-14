% Thi script analyses the performance of a reduced Walker constellation
% under low satellite budget requirements.

% This script replicates results shown in Figures 27-28 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

%% Parameters

Constants = initialiseAstronomicalConstants();
Parameters.nOrb = 15;
Parameters.nSatOrb = 27;
Parameters.h = 460;
Parameters.in = deg2rad(77.1);
Parameters.raan = 0;
Parameters.nSatOrbActual = 5;
areOrbitsOverlapped = 0;
Parameters.dtheta = deg2rad(8.5) - ...
    areOrbitsOverlapped * 2*pi * Parameters.nSatOrbActual / Parameters.nSatOrb;
sizeThetaVec = 100;
dirSun = [-1 0 0]';
T = computeOrbitalPeriod(Constants.R_E + Parameters.h, Constants.MU_E);
%deltaT = T / (Parameters.nSatOrb);
deltaT = T / 100;
tspan = 0:deltaT:T;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e3; % substitute 1e4
CoverageOptions.coverageCount = 3; % Triple coverage
R = Constants.R_E + 100;
deltaR = 2000;
dispFlag = 0;

flagOverlappedAnalysis = 1; % 1 if you want to test effect of satellite overlapping
flagParametricAnalysis = 1; % 1 if you want a parametric analysis on the number of satellites in the constellation

areOrbitsOverlappedVec = [1, 0];
legendAreOrbitsOverlapped = ["overlapped", "not overlapped"];

nSatOrbActualVec = [3, 5, 7, 9]; % Overlapped case
for k = 1:size(nSatOrbActualVec,2)
    legendNSatOrbActual(k) = "N_{sat/o} = "+string(nSatOrbActualVec(k));
end

%% Analysis - overlapping v not overlapping

if flagOverlappedAnalysis
    instantaneousCoverageMat = [];
    cumulativeCoverageMat = [];
    tic
    for k = 1:size(areOrbitsOverlappedVec, 2)
        Parameters.dtheta = deg2rad(8.5) - ...
        areOrbitsOverlappedVec(k) * 2*pi * Parameters.nSatOrbActual / Parameters.nSatOrb;
        [instantaneousCoverageVec, cumulativeCoverageVec] = ...
            coverageAnalysis(Parameters, tspan, SampleOptions, dirSun, R, deltaR, ...
            CoverageOptions, Constants, dispFlag);
        instantaneousCoverageMat(:,k) = instantaneousCoverageVec';
        cumulativeCoverageMat(:,k) = cumulativeCoverageVec';
        disp(k+" / "+size(areOrbitsOverlappedVec, 2)+", time: "+toc)
    end
    
    figure()
    
    subplot(2,1,1)
    hold on
    for k = 1:size(areOrbitsOverlappedVec, 2)
        plot(tspan/T, 100*instantaneousCoverageMat(:,k)', "LineWidth", 2)
    end
    %lgd = legend(legendAreOrbitsOverlapped);
    fontsize(gca, 12, 'points')
    ylabel("inst. coverage [%]", "FontSize", 12)
    
    subplot(2,1,2)
    hold on
    for k = 1:size(areOrbitsOverlappedVec, 2)
        plot(tspan/T, 100*cumulativeCoverageMat(:,k)', "LineWidth", 2)
    end
    lgd = legend(legendAreOrbitsOverlapped);
    lgd.Location = "southeast";
    fontsize(gca, 12, 'points')
    ylabel("cumulative coverage [%]", "FontSize", 12)
    xlabel("time [T]", "FontSize", 12)
end

%% Analysis - number of satellites per orbit

if flagParametricAnalysis
    instantaneousCoverageMat = [];
    cumulativeCoverageMat = [];
    tic
    for k = 1:size(nSatOrbActualVec, 2)
        Parameters.nSatOrbActual = nSatOrbActualVec(k);
        Parameters.dtheta = deg2rad(8.5) - ...
        2*pi * Parameters.nSatOrbActual / Parameters.nSatOrb;
        [instantaneousCoverageVec, cumulativeCoverageVec] = ...
            coverageAnalysis(Parameters, tspan, SampleOptions, dirSun, R, deltaR, ...
            CoverageOptions, Constants, dispFlag);
        instantaneousCoverageMat(:,k) = instantaneousCoverageVec';
        cumulativeCoverageMat(:,k) = cumulativeCoverageVec';
        disp(k+" / "+size(nSatOrbActualVec, 2)+", time: "+toc)
    end
    
    figure()
    
    subplot(2,1,1)
    hold on
    for k = 1:size(nSatOrbActualVec, 2)
        plot(tspan/T, 100*instantaneousCoverageMat(:,k)', "LineWidth", 2)
    end
    %lgd = legend(legendNSatOrbActual);
    fontsize(gca, 12, 'points')
    ylabel("inst. coverage [%]", "FontSize", 12)
    
    subplot(2,1,2)
    hold on
    for k = 1:size(nSatOrbActualVec, 2)
        plot(tspan/T, 100*cumulativeCoverageMat(:,k)', "LineWidth", 2)
    end
    lgd = legend(legendNSatOrbActual);
    lgd.Location = "southeast";
    fontsize(gca, 12, 'points')
    ylabel("cumulative coverage [%]", "FontSize", 12)
    xlabel("time [T]", "FontSize", 12)
end

%% Plots

% figure()
% subplot(2,1,1)
% plot(tspan/T, 100*instantaneousCoverageVec)
% subplot(2,1,2)
% plot(tspan/T, 100*cumulativeCoverageVec)

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
function flag = isVisible(rObs, rTrg, dirSun, R, deltaR)
    flag    = ~( ( ~isWithinRange(rObs, rTrg, deltaR) ) || ... 
            isScatteringWeak(rObs, rTrg, dirSun) || isEclipsed(rTrg, dirSun, R) || ...
            ~isAboveTheHorizon(rObs, rTrg, R) );
end

% Coverage analysis function 
function [instantaneousCoverageVec, cumulativeCoverageVec] = coverageAnalysis(Parameters, tspan, SampleOptions, dirSun, R, deltaR, CoverageOptions, Constants, dispFlag)
    % Initialise constellation
    Constellation = initialiseLowBudgetConstellation(Parameters, Constants);
    %nSat = size(Constellation, 2);
    
    % Propagate constellation
    nT = size(tspan, 2);
    Constellation = propagateConstellation(tspan, Constellation, Constants.MU_E);
    
    % Initialise Monte Carlo sample (fixed for the simulation)
    rTrgMat = initialiseSampleShell(SampleOptions);
    
    % Compute istantaneous coverage at each time instant
    instantaneousCoverageVec = nan(1, nT);
    if dispFlag
        tic
    end
    for kT = 1:nT
        rObsMat = convertStructToMat(Constellation, kT);
        instantaneousCoverageVec(kT) = computeInstantaneousCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
        if dispFlag
            disp("Instantaneous coverage computation, "+kT+" / "+nT+", time: "+toc)
        end
    end
    
    % Compute cumulative coverage
    cumulativeCoverageVec = computeCumulativeCoverageVec(Constellation, nT, rTrgMat, dirSun, R, deltaR, CoverageOptions, dispFlag);
end

% Cumulative coverage function
function cumulativeCoverageVec = computeCumulativeCoverageVec(Constellation, nT, rTrgMat, dirSun, R, deltaR, CoverageOptions, displayFlag)
    cumulativeCoverageVec = nan(1, nT);
    nTrg = size(rTrgMat, 1);
    hasBeenTripleVisible = zeros(1, nTrg);
    if displayFlag
        tic
    end
    for kT = 1:nT
        rObsMat = convertStructToMat(Constellation, kT);
        nObs = size(rObsMat, 1);
        coverage = 0;
        for iTrg = 1:nTrg
            rTrg = rTrgMat(iTrg,:)';
            count = 0;
            iObs = 1;
            while count < CoverageOptions.coverageCount && iObs < nObs
                rObs = rObsMat(iObs,:)';
                count = count + isVisible(rObs, rTrg, dirSun, R, deltaR);
                iObs = iObs + 1;
            end
            isTripleVisibile = (count == CoverageOptions.coverageCount);
            hasBeenTripleVisible(iTrg) = (hasBeenTripleVisible(iTrg) || isTripleVisibile);
            coverage = coverage + hasBeenTripleVisible(iTrg);
        end
        cumulativeCoverageVec(kT) = coverage / nTrg;
        if displayFlag
            disp("Timestep: "+kT+" / "+nT+", Time: "+toc)
        end
    end
end

% Define instantaneous coverage function
function coverage = computeInstantaneousCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions)
    nObs = size(rObsMat, 1);
    nTrg = size(rTrgMat, 1);
    coverage = 0;
    for iTrg = 1:nTrg
        rTrg = rTrgMat(iTrg,:)';
        count = 0;
        iObs = 1;
        while count < CoverageOptions.coverageCount && iObs < nObs
            rObs = rObsMat(iObs,:)';
            count = count + isVisible(rObs, rTrg, dirSun, R, deltaR);
            iObs = iObs + 1;
        end
        coverage = coverage + (count == CoverageOptions.coverageCount);
    end
    coverage = coverage / nTrg;
end

% Define average coverage function
function coverage = computeCoverage(timeVec, Constellation, rTrgMat, R, deltaR, CoverageOptions, Constants)
    Constellation = propagateConstellation(timeVec, Constellation, Constants);
    totalCoverage = 0;
    for timeStep = 1:size(timeVec, 2)
        rObsMat = convertStructToMat(Constellation, timeStep);
        totalCoverage = totalCoverage + computeInstantaneousCoverage(rObsMat, rTrgMat, R, deltaR, CoverageOptions);
    end
    coverage = totalCoverage / size(timeVec, 2);
end