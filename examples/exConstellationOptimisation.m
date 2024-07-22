% This example performs a global optimisation task to find the
% constellation parameters maximising coverage, given a budget in terms of
% number of constellation satellites.

% Ex. nSc = 400 -> 8.0000   50.0000    1.2428    0.0625    0.0951, 63.5%
% 84.0000    5.0000    1.1775    0.0199    0.3320, 59.5%
% 8.0000   49.0000    1.3312    0.3694    0.0213, 64.95%

clc
clear
close all

%% Settings

Constants = initialiseAstronomicalConstants();
dirSun = [-1, 0, 0]';
R = Constants.R_E + 100;
deltaR = 2000;
h = 460;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e4;
CoverageOptions.coverageCount = 3;
OptimisationOptions     = optimoptions('ga', 'PopulationSize', 40, 'PlotFcn', @gaplotbestf,...
                        'MaxStallGenerations', 40, 'MaxGenerations', 60, "FunctionTolerance", 5*1e-3);
nSc = 400;
deltaNSc = 10;
lowerBound = [3, floor(sqrt(nSc - deltaNSc)), deg2rad(30), 0, 0]; % nOrb, nSatOrb, in, raan, dtheta
upperBound = [ceil(sqrt(nSc + deltaNSc)), nSc + deltaNSc, deg2rad(80), 2*pi, 2*pi]; % nOrb, nSatOrb, in, raan, dtheta
integerVar = [1 2]; % integer parameters (nOrb, nSatOrb)
nVars = 5;
nonCon = @(p) nonlinearConstraints(p, nSc, deltaNSc);
SurrogateOptimisationOptions    = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
                                'MaxFunctionEvaluations', 120);

%% Optimisation

fitnessFunction = @(p) fitnessResidualCoverage(p, dirSun, R, deltaR, h, Constants, SampleOptions, CoverageOptions);
[optimalParameters, optimalResidualCoverage]    = ga(fitnessFunction, ...
        nVars, [], [], [], [], lowerBound, upperBound, nonCon, integerVar, ...
        OptimisationOptions);

% surrogateFunctionHandle = @(p) surrogateFunction(p, fitnessFunction, nonCon);
% [optimalParameters, optimalResidualCoverage]    = surrogateopt(surrogateFunctionHandle, ...
%          lowerBound, upperBound, integerVar, ...
%          SurrogateOptimisationOptions);

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
            isBlinded(rObs, rTrg, dirSun) || isEclipsed(rTrg, dirSun, R) || ...
            isObstructed(rObs, rTrg, R) );
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

% Define optimisation function to be minimised
function residualCoverage   = fitnessResidualCoverage(p, dirSun, R, deltaR, h, ...
                            Constants, SampleOptions, CoverageOptions)
    % Extract parameters
    Parameters.nOrb = p(1);
    Parameters.nSatOrb = p(2);
    Parameters.in = p(3);
    Parameters.raan = p(4);
    Parameters.dtheta = p(5);
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

% Define nonlinear constraints
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

function f = surrogateFunction(p, fitnessFunction, nonCon)
    f.Fval = fitnessFunction(p);
    f.Ineq = nonCon(p);
end