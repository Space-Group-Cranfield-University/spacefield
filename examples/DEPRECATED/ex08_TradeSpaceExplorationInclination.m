% This script examines how coverage varies with the inclination of the
% constellation at fixed parameters.
%
% This script replicates results shown in Fig. 18 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
close all
clc

Constants = initialiseAstronomicalConstants();

%% Settings

nConstellations = 40; % Number of initialised constellations (500 <-> 15'). Substitute 200
nMc = 1e3; % Number of Monte Carlo points. Substitute 1e4
Parameters.h = 460; % Constellation altitude
%nInterval = [100 1000];
nOrbVec = [6 8 10];
Parameters.nSatOrb = 50;
Parameters.raan = 0;
Parameters.dtheta = deg2rad(5);
inVec = deg2rad(20:2:160); % Interval of allowed inclinations
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = nMc;
dirSun = [-1 0 0]';
R = Constants.R_E + 100;
deltaR = 2000;
CoverageOptions.coverageCount = 3;

%% Trade space exploration

%inVec = nan(1, nConstellations);
coverageMat = nan(size(nOrbVec, 2), size(inVec, 2));
tic
for k = 1:size(nOrbVec, 2)
    Parameters.nOrb = nOrbVec(k);
    for iConstellations = 1:size(inVec, 2)
        % Initialise constellation
        Parameters.in = inVec(iConstellations);
        Constellation = initialiseConstellation(Parameters, Constants);
        rObsMat = convertStructToMat(Constellation);
    
        % Initialise Monte Carlo points
        rTrgMat = initialiseSampleShell(SampleOptions);
    
        % Compute coverage
        coverage = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
        %residualCoverage = 1 - coverage;
    
        % Allocate results
        %inVec(iConstellations) = Parameters.in;
        coverageMat(k, iConstellations) = coverage;
    
        if mod(iConstellations, 20) == 0
            disp("Constellations tested: "+iConstellations+" / "+size(inVec, 2)+". Time: "+toc)
        end
    end
end

%[inPareto, coveragePareto] = estimateParetoFrontier(inVec, coverageVec);

%% Plots

figure()
% plot(rad2deg(inVec), coverageMat(1,:)*100, "LineWidth", 1)
scatter(rad2deg(inVec), coverageMat(1,:)*100, "v", "filled", "LineWidth", 1)
hold on
scatter(rad2deg(inVec), coverageMat(2,:)*100, "filled", "LineWidth", 1)
scatter(rad2deg(inVec), coverageMat(3,:)*100, "^", "filled", "LineWidth", 1)
%plot(rad2deg(inPareto), coveragePareto*100, "LineWidth", 2)
%lgd = legend(["Local points","Best performance"]);
%lgd.Location = "northwest";
legend(["N = 300","N = 400","N = 500"])
fontsize(gca, 12, 'points')
xlabel("constellation inclination [deg]", "FontSize", 14)
ylabel("triple coverage [%]", "FontSize", 14)
title("Triple coverage against Inclination", "FontSize", 12)

%% Functions

function [inPareto, coveragePareto] = estimateParetoFrontier(inVec, coverageVec)
    [inVecSorted, sortIndex] = sort(inVec);
    coverageVecSorted = coverageVec(sortIndex);
    coveragePareto = coverageVecSorted(1);
    inPareto = inVecSorted(1);
    for k = 2:length(coverageVecSorted)
        if inVecSorted(k) == inPareto(end) % Case where two coverage values have the same number of satellites
            if coverageVecSorted(k) > coveragePareto(end)
                coveragePareto(end) = coverageVecSorted(k);
                inPareto(end) = inVecSorted(k);
            end
        elseif coverageVecSorted(k) > coveragePareto(end)
            coveragePareto = [coveragePareto, coverageVecSorted(k)];
            inPareto = [inPareto, inVecSorted(k)];
        end
    end
end

function Parameters = initialiseSampleParameters(nInterval, inInterval, fixedParameters)

    % Function sampling the parameters of a constellation. For real-valued
    % parameters, the function samples uniformly over the given intervals
    % or, in the case of altitude, sets it to the given fixed altitude
    % parameter. For the integer-valued parameters it works as follows:
    %
    % Input:    budget interval. [100 1000] means we have a budget
    %           between 100 and 1000 satellites.
    %
    %   1)  extracts a random budget N from the budget interval using randi()
    %   2)  factors N using factor()
    %   3)  if N is prime, set nSatOrb = N, nOrb = 1.
    %   4)  else, select a random and randomly sized subset of the factors
    %       of N. If the product of the subset is greater than sqrt(N),
    %       assign it to nSatOrb,and N/nSatOrb to nOrb. Else, assign it to
    %       nOrb and N/nOrb to nSatOrb.
    %
    % Assumption: nSatOrb >= nOrb   
    
    Parameters.nSat = randi(nInterval);
    sqrtN = sqrt(Parameters.nSat);
    factors = factor(Parameters.nSat);
    if size(factors, 2) > 1
        nSample = randi(size(factors, 2));
        nOrbFactors = randsample(factors, nSample);
        factorProduct = prod(nOrbFactors);
        if factorProduct >= sqrtN
            Parameters.nSatOrb = factorProduct;
            Parameters.nOrb = Parameters.nSat / Parameters.nSatOrb;
        else
            Parameters.nOrb = factorProduct;
            Parameters.nSatOrb = Parameters.nSat / Parameters.nOrb;
        end
    else % Case: prime number
        Parameters.nOrb = 1;
        Parameters.nSatOrb = Parameters.nSat / Parameters.nOrb;
    end

    Parameters.in = inInterval(1) + rand * (inInterval(2) - inInterval(1));
    Parameters.h = fixedParameters.h;
    Parameters.raan = rand * 2 * pi / Parameters.nOrb;
    Parameters.dtheta = rand * 2 * pi / Parameters.nSatOrb;
end

% Transform constellation structure in matrix
function rObsMat = convertStructToMat(Constellation)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x0(1:3,1);
    end
end

% Define visibility function
function flag = isVisible(rObs, rTrg, dirSun, R, deltaR)
    flag = isEclipsed(rTrg, dirSun, R);
    if flag == 0
        flag = ~isWithinRange(rObs, rTrg, deltaR);
        if flag == 0
            flag = ~isAboveTheHorizon(rObs, rTrg, R);
            if flag == 0
                flag = isScatteringWeak(rObs, rTrg, dirSun);
            end
        end
    end
    flag = ~flag;
end

% Define instantaneous coverage function
function coverage = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions)
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