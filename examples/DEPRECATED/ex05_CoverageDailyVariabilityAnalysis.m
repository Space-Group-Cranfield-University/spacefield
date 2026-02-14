% This script analyses how the coverage varies daily due to the motion of
% the constellation itself.

% This script replicates results shown in Figures 13-14 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

%% Settings

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
SampleOptions.N = 1e3; % substitute 1e4
dirSun = [-1 0 0]';
nSatOrbVec = 10:5:50;
nTime_1 = 50;
nTime_2 = 20;
test_1 = 1; % 1 if you want to test time variability for test constellation
test_2 = 1; % 1 if you want to test how time variability varies with nSatOrb

%% Processing

orbitalPeriod = computeOrbitalPeriod(Constants.R_E + ConstellationParameters.h, Constants.MU_E);
timeVec_1 = linspace(0, orbitalPeriod, nTime_1);
timeVec_2 = linspace(0, orbitalPeriod, nTime_2);

%% Test 1 (time variability for test constellation)

if test_1
    % Initialise constellation
    Constellation = initialiseConstellation(ConstellationParameters, Constants);
    
    % Propagate constellation
    Constellation = propagateConstellation(timeVec_1, Constellation, Constants.MU_E);
    
    % Initialise struct
    coverageVec = nan(1, nTime_1);
    tic
    for kTime = 1:nTime_1
        % Initialise Monte Carlo sample
        rTrgMat = initialiseSampleShell(SampleOptions);

        rObsMat = convertStructToMat(Constellation, kTime);
        coverageVec(kTime)  = computeInstantaneousCoverage(rObsMat, rTrgMat,...
                            dirSun, R, deltaR, CoverageOptions);
        disp(kTime+" / "+nTime_1+", time: "+toc)
    end

    % Plot
    figure()
    plot(timeVec_1./orbitalPeriod, coverageVec*100, "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("time [T]", "FontSize", 14)
    ylabel("triple coverage [%]", "FontSize", 14)
    title("Triple coverage against Time", "FontSize", 12)
end

%% Test 2 (parametric analysis time variability vs test constellation)

if test_2
    
    coverageMat = nan(size(nSatOrbVec, 2), nTime_2);
    deltaCoverageVec = nan(1, size(nSatOrbVec, 2));
    tic
    for kNSatOrb = 1:size(nSatOrbVec, 2)
        ConstellationParameters.nSatOrb = nSatOrbVec(kNSatOrb);
        coverageVec = computeTimeCoverageWithoutResampling(timeVec_2, dirSun, R, deltaR,...
                    CoverageOptions, SampleOptions, ConstellationParameters, Constants);
        deltaCoverageVec(kNSatOrb) = max(coverageVec) - min(coverageVec);
        coverageMat(kNSatOrb, :) = coverageVec;
        disp("nSatOrb: "+ConstellationParameters.nSatOrb+" , time: "+toc)
    end
    meanCoverageVec = mean(coverageMat, 2);

    % Plot
    figure()
    plot(nSatOrbVec, deltaCoverageVec*100, "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("number of satellites per orbit", "FontSize", 14)
    ylabel("max cov_3 - min cov_3 [%]", "FontSize", 14)
    title("Delta triple coverage against Number of satellites per orbit", "FontSize", 12)

    figure()
    semilogy(nSatOrbVec, 100*deltaCoverageVec./meanCoverageVec', "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("number of satellites per orbit", "FontSize", 14)
    ylabel("\Deltacov_3 / cov_{3} [%]", "FontSize", 14)
    title("Delta triple coverage over Average triple coverage", "FontSize", 12)

    figure()

    subplot(2,1,1)
    semilogy(nSatOrbVec, 100*meanCoverageVec', "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("number of satellites per orbit", "FontSize", 14)
    ylabel("cov_{3,avg} [%]", "FontSize", 14)
    title("Average triple coverage", "FontSize", 12)

    subplot(2,1,2)
    semilogy(nSatOrbVec, 100*deltaCoverageVec./meanCoverageVec', "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("number of satellites per orbit", "FontSize", 14)
    ylabel("\Deltacov_3 / cov_{3,avg} [%]", "FontSize", 14)
    title("Delta triple coverage over Average triple coverage", "FontSize", 12)
end

%% Functions

% Transform constellation structure in position matrix
function rObsMat = convertStructToMat(Constellation, kTime)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x(kTime, 1:3)';
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
function coverage = computeInstantaneousCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions)
    nObs = size(rObsMat, 1);
    integrandFunction = @(r) isVisibleToConstellation(r, rObsMat, nObs, dirSun, R, deltaR, CoverageOptions);
    coverage = integrateMonteCarlo(rTrgMat, integrandFunction);
end

% Single test function
function coverageVec = computeTimeCoverageWithoutResampling(timeVec, dirSun, R, deltaR, CoverageOptions, SampleOptions, ConstellationParameters, Constants)
    nTime = size(timeVec, 2);

    % Initialise Monte Carlo sample
    rTrgMat = initialiseSampleShell(SampleOptions);

    % Initialise constellation
    Constellation = initialiseConstellation(ConstellationParameters, Constants);
    
    % Propagate constellation
    Constellation = propagateConstellation(timeVec, Constellation, Constants.MU_E);
    
    % Initialise struct
    coverageVec = nan(1, nTime);
    tic
    for kTime = 1:nTime
        rObsMat = convertStructToMat(Constellation, kTime);
        coverageVec(kTime) = computeInstantaneousCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
        disp(kTime+" / "+nTime+", time: "+toc)
    end
end