% This script examines how coverage varies with the number of constellation
% satellites with random parameters. The results are presented in terms of
% residual coverage, i.e. 1 - coverage, to show the Pareto front.

clear
close all

Constants = initialiseAstronomicalConstants();

%% Settings

nConstellations = 500; % Number of initialised constellations (500 <-> 15')
nMc = 1e4; % Number of Monte Carlo points
fixedParameters.h = 460; % Constellation altitude
nInterval = [100 600]; % Interval of possible number of constellation satellites
inInterval = deg2rad([20 80]); % Interval of allowed inclinations
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = nMc;
dirSun = [-1 0 0]';
R = Constants.R_E + 100;
deltaR = 2000;
CoverageOptions.coverageCount = 3;

%% Trade space exploration

nSatVec = nan(1, nConstellations);
residualCoverageVec = nan(1, nConstellations);
tic
for iConstellations = 1:nConstellations

    % Initialise constellation
    Parameters = initialiseSampleParameters(nInterval, inInterval, fixedParameters);
    Constellation = initialiseConstellation(Parameters, Constants);
    rObsMat = convertStructToMat(Constellation);

    % Initialise Monte Carlo points
    rTrgMat = initialiseSampleShell(SampleOptions);

    % Compute coverage
    coverage = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
    residualCoverage = 1 - coverage;

    % Allocate results
    nSatVec(iConstellations) = Parameters.nSat;
    residualCoverageVec(iConstellations) = residualCoverage;

    if mod(iConstellations, 20) == 0
        disp(iConstellations)
        disp(toc)
    end

end

%% Plots

figure()
scatter(nSatVec, residualCoverageVec*100)
xlabel("Number of constellation satellites")
ylabel("Residual triple coverage [%]")
title("Residual coverage VS Number of constellation satellites")

%% Functions

function Parameters = initialiseSampleParameters(nInterval, inInterval, fixedParameters)

    Parameters.nSat = randi(nInterval);
    factors = factor(Parameters.nSat);
    if size(factors, 2) > 1
        factorsTemp = factors;
        factors(factors == 2) = [];
        if size(factors, 2) == 0
            Parameters.nSatOrb = prod(factorsTemp);
            Parameters.nOrb = Parameters.nSat / Parameters.nSatOrb;
        elseif size(factors, 2) > 1
            nSample = randi(size(factors, 2));
            nOrbFactors = randsample(factors, nSample);
            Parameters.nOrb = prod(nOrbFactors);
            Parameters.nSatOrb = Parameters.nSat / Parameters.nOrb;
        else
            Parameters.nOrb = factors;
            Parameters.nSatOrb = Parameters.nSat / Parameters.nOrb;     
        end
    elseif randi([0, 1])
        Parameters.nOrb = factors;
        Parameters.nSatOrb = Parameters.nSat / Parameters.nOrb;
    else
        Parameters.nSatOrb = factors;
        Parameters.nOrb = Parameters.nSat / Parameters.nSatOrb;
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
            flag = isObstructed(rObs, rTrg, R);
            if flag == 0
                flag = isBlinded(rObs, rTrg, dirSun);
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