% This example performs a global optimisation task to find the
% constellation parameters maximising coverage, given a budget in terms of
% number of constellation satellites.

% This script replicates results shown in Fig. 19 at https://doi.org/10.1016/j.actaastro.2025.02.019

clc
clear
close all

%% Settings

Constants = initialiseAstronomicalConstants();
dirSun = [-1, 0, 0]';
R = Constants.R_E + 100;
deltaR = 2000;
h = 460;
nSc = 400; % Budget
deltaNSc = 5; % Acceptable variation from budget
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e3; % substitute 1e4
CoverageOptions.coverageCount = 3;

% Random Search options
nCalls = 10; % substitute 100
nRuns = 10;
plotRandomSearch = 1;
Bounds.nSat = [nSc - deltaNSc, nSc + deltaNSc];
Bounds.in = deg2rad([75, 80]);
Bounds.h = h;
dispOptions = 1;
plotOptions = 1;

% GA options
OptimisationOptions     = optimoptions('ga', 'PopulationSize', 40, 'PlotFcn', @gaplotbestf,...
                        'MaxStallGenerations', 40, 'MaxGenerations', 60, "FunctionTolerance", 5*1e-3);
lowerBound = [1, floor(sqrt(nSc - deltaNSc)), deg2rad(75), 0]; % nOrb, nSatOrb, in, raan, dtheta
upperBound = [ceil(sqrt(nSc + deltaNSc)), nSc + deltaNSc, deg2rad(80), 2*pi]; % nOrb, nSatOrb, in, raan, dtheta
integerVar = [1 2]; % integer parameters (nOrb, nSatOrb)
nVars = 4;
nonCon = @(p) nonlinearConstraints(p, nSc, deltaNSc);

randomSearch = 1;
geneticAlgorithm = 0;

%% Optimisation

if randomSearch
    bestCoverageMat = nan(nRuns, nCalls);
    tic
    for k = 1:nRuns
        costFunction = @(p) costFunctionCoverage(p, dirSun, R, deltaR, h, ...
                                Constants, SampleOptions, CoverageOptions);
        [bestCoverage, bestParameters, bestCoverageVec, avgCoverageVec] = ...
            executeRandomSearch(costFunction, Bounds, nCalls, dispOptions);
        bestCoverageMat(k, :) = bestCoverageVec;
        if dispOptions
            disp("Run: "+k+" / "+nRuns+", time: "+toc)
        end
    end
    avgCoverageMat = mean(bestCoverageMat, 1);
    if plotOptions
        if nRuns == 1
            figure()
            hold on
            plot(1:nCalls, bestCoverageVec*100, "LineWidth", 2)
            plot(1:nCalls, avgCoverageVec*100, "LineWidth", 2)
            lgd = legend(["best", "mean"]);
            lgd.Location = "southeast";
            fontsize(gca, 12, 'points')
            xlabel("random search calls", "FontSize", 14)
            ylabel("triple coverage [%]", "FontSize", 14)
            title("Random search of best triple coverage", "FontSize", 12)
        else
            figure()
            hold on

            for k = 1:nRuns
                if k == 1
                    p_1 = stairs(1:nCalls, bestCoverageMat(k,:) * 100, "c");
                else
                    stairs(1:nCalls, bestCoverageMat(k,:) * 100, "c");
                end
            end
            p_avg = plot(1:nCalls, avgCoverageMat*100, "r", "LineWidth", 2);
            lgd = legend([p_1,p_avg], "search run", "mean");
            lgd.Location = "southeast";
            fontsize(gca, 12, 'points')
            xlabel("random search calls", "FontSize", 14)
            ylabel("best triple coverage [%]", "FontSize", 14)
            title("Random search of best triple coverage", "FontSize", 12)
        end
    end
end

if geneticAlgorithm
    fitnessFunction = @(p) fitnessResidualCoverage(p, dirSun, R, deltaR, h,...
                    Constants, SampleOptions, CoverageOptions);
    [optimalParameters, optimalResidualCoverage]    = ga(fitnessFunction, ...
                                                    nVars, [], [], [], [], ...
                                                    lowerBound, upperBound, ...
                                                    nonCon, integerVar, ...
                                                    OptimisationOptions);
end

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
function flag = isVisible(rObs, rTrg, dirSun, R, deltaR)
    flag    = ~( ( ~isWithinRange(rObs, rTrg, deltaR) ) || ... 
            isScatteringWeak(rObs, rTrg, dirSun) || isEclipsed(rTrg, dirSun, R) || ...
            ~isAboveTheHorizon(rObs, rTrg, R) );
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

% Define optimisation function to be maximised
function coverage   = costFunctionCoverage(Parameters, dirSun, R, deltaR, h, ...
                            Constants, SampleOptions, CoverageOptions)
    % Modify parameters
    Parameters.raan = 0;

    % Initialise constellation
    Constellation = initialiseConstellation(Parameters, Constants);
    rObsMat = convertStructToMat(Constellation);

    % Initialise Monte Carlo points
    rTrgMat = initialiseSampleShell(SampleOptions);

    % Compute coverage
    coverage = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
end

% Random search function
function [bestCoverage, bestParameters, bestCoverageVec, avgCoverageVec] = ...
    executeRandomSearch(costFunction, Bounds, nCalls, dispOptions)
    
    avgCoverageVec = nan(1, nCalls);
    bestCoverageVec = nan(1, nCalls);
    bestCoverage = 0;
    bestParameters = nan;
    tic
    for k = 1:nCalls
        parameters = initialiseSampleParameters(Bounds.nSat, Bounds.in, Bounds);
        coverage = costFunction(parameters);
        if coverage > bestCoverage
            bestCoverage = coverage;
            bestParameters = parameters;
        end
        bestCoverageVec(k) = bestCoverage;
        if k > 1
            avgCoverageVec(k) = ((k-1) * avgCoverageVec(k-1) + coverage)/k;
        else
            avgCoverageVec(k) = coverage;
        end
        if dispOptions
            disp("Call: "+k+" / "+nCalls+", time: "+toc)
        end
    end

end

% Define optimisation function to be minimised
function residualCoverage   = fitnessResidualCoverage(p, dirSun, R, deltaR, h, ...
                            Constants, SampleOptions, CoverageOptions)
    % Extract parameters
    Parameters.nOrb = p(1);
    Parameters.nSatOrb = p(2);
    Parameters.in = p(3);
    Parameters.raan = 0;
    Parameters.dtheta = p(4);
    Parameters.h = h;

    % Initialise constellation
    Constellation = initialiseConstellation(Parameters, Constants);
    rObsMat = convertStructToMat(Constellation);

    % Initialise Monte Carlo points
    rTrgMat = initialiseSampleShell(SampleOptions);

    % Compute coverage
    coverage = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
    residualCoverage = 1 - coverage;
end

% Define nonlinear constraints for GA
function [c, ceq] = nonlinearConstraints(p, nSc, deltaNSc)
    % p contains nOrb, nSatOrb, in, raan, dtheta
    % equality constraint expressed in terms of tolerance
    tol = deltaNSc;
    ceq = [];
    c   = [ nSc - p(1) * p(2) - tol; 
            p(1) * p(2) - nSc - tol;
            p(4) - 2 * pi / p(1); 
            p(5) - 2 * pi / p(2)    ];
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