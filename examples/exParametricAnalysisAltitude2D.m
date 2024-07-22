% This script provides a parametric analysis for the average coverage of a
% two-dimensional constellation. The parameter is the constellation
% altitude. The coverage is computed using Monte Carlo integration. The
% constellation is assumed to be a set of equispaced satellites on a
% circular equatorial orbit. We consider the triple coverage, i.e. the
% average percentage of volume that is covered by at least three observers.

clear
close all

%% Settings

% Set constants
H_LEO_LOW = 400; % Lower altitude of LEO shell
H_LEO_HIGH = 1800; % Higher altitude of LEO shell
Constants = initialiseAstronomicalConstants();
R = Constants.R_E_NORAD + 100; % Radius of the Earth plus atmospheric height
deltaR = 2000; % Visibility range

% Set simulation parameters
nMC = 1e4; % Number of Monte Carlo points
nH = 20; % Size of altitude parameter vector
stepN = 10; % Size of number of constellation satellites vector
SampleOptions.dim = "2D";
CoverageOptions.coverageCount = 3; % We consider triple coverage
deltaT = 100; % 30 seconds timestep

% Set constellation parameters
hVec = linspace(H_LEO_LOW, H_LEO_HIGH, nH);
nSatVec = 50:stepN:90;

% Set Sun direction
dirSun = [-1 0]';

%% Process settings
SampleOptions.R_1 = Constants.R_E + H_LEO_LOW;
SampleOptions.R_2 = Constants.R_E + H_LEO_HIGH;
SampleOptions.N = nMC;
nN = size(nSatVec, 2);
for iN = 1:nN
    legendVect(iN) = "nSat: "+string(nSatVec(iN));
end

%% Sample target points
rTrgMat = initialiseSampleShell(SampleOptions);

%% Parametric analysis
coverageResults = nan(nN, nH);
for iN = 1:nN
    nSatConstellation = nSatVec(1,iN);
    for iH = 1:nH
        constellationAltitude = hVec(1,iH);
        coverage    = computeAverageCoverage(deltaT, constellationAltitude,...
                    nSatConstellation, rTrgMat, dirSun, R, deltaR, CoverageOptions, ...
                    Constants);
        coverageResults(iN, iH) = coverage;
        disp("Done altitude: "+iH)
    end
    disp("Done nSat: "+iN)
end

%% Plots
figure()
hold on
for iN = 1:nN
    plot(hVec, coverageResults(iN,:)*100)
end
xlabel("Altitude [km]")
ylabel("Triple coverage [%]")
title("Triple coverage VS Constellation altitude (2D case)")
legend(legendVect)

%% Coverage functions

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
function coverage   = computeAverageCoverage(deltaT, constellationAltitude,...
                    nSatConstellation, rTrgMat, dirSun, R, deltaR, CoverageOptions, ...
                    Constants)
    t = 0; % starting time
    nT = 0; % number of timesteps

    % Compute period of the constellation orbit
    periodConstellation = computeOrbitalPeriod(Constants.R_E + ...
                        constellationAltitude, Constants.MU_E);
    
    % Compute final time of propagation (the motion of the constellation is
    % periodic afterwards)
    tFinal = periodConstellation / nSatConstellation;
    
    rObsMat = zeros(nSatConstellation, 2);
    rConstellation = Constants.R_E + constellationAltitude;
    coverage = 0;
    while t < tFinal
        % Propagate constellation
        for iObs = 1:nSatConstellation
            phase = (iObs-1) * 2*pi / nSatConstellation;
            rObsMat(iObs, 1) = rConstellation * cos(2*pi*t/periodConstellation + phase);
            rObsMat(iObs, 2) = rConstellation * sin(2*pi*t/periodConstellation + phase);
        end

        coverage = coverage + computeInstantaneousCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
        t = t + deltaT;
        nT = nT + 1;
    end

    if nT > 0
        coverage = coverage / nT;
    else
        coverage = -1; % manages division by zero
    end
end