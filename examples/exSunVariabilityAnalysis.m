% This script analyses how the yearly Sun position affects the coverage of
% a given constellation, both for fixed and variable constellation orbit
% inclination.

clear
clc
close all

Constants = initialiseAstronomicalConstants;
dirSun0 = [-1 0 0]';
%dirSun1 = [0 -cosd(Constants.ECL) sind(Constants.ECL)]';
nTime = 13;
timeVec = linspace(0, 52*7*86400, nTime);
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
SampleOptions.N = 5*1e4; % 1e4 too little to clearly see dependency on Sun
inVec = 20:10:80; % vector of inclinations
inVec = deg2rad(inVec);
nIn = size(inVec,2);
yearlyVariabilityAnalysis = 0;
parametricAnalysis = 1;

%% Sun propagation

dirSunMat = propagateDirSun(timeVec, dirSun0, Constants);

%% Compute coverage

% Initialise Monte Carlo sample
rTrgMat = initialiseSampleShell(SampleOptions);

% Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);
rObsMat = convertStructToMat(Constellation);

% Coverage computation
if yearlyVariabilityAnalysis
    coverageVec = nan(1, nTime);
    tic
    for iTime = 1:nTime
        dirSun = dirSunMat(iTime, :)';
        coverageVec(iTime) = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
        disp("Timestep: "+iTime+" / "+nTime+", time: "+toc)
    end  
    figure()
    plot(timeVec/86400, coverageVec*100)
    title("Triple coverage throughout the year")
    xlabel("time [days]")
    ylabel("coverage [%]")
end

% Parametric analysis on constellation inclination

skip = 1;
if parametricAnalysis
    tic
    coverageMat = nan(nIn, nTime / skip);
    for iIn = 1:nIn
        ConstellationParameters.in = inVec(iIn);
        Constellation = initialiseConstellation(ConstellationParameters, Constants);
        rObsMat = convertStructToMat(Constellation);
        for iTime = 1:skip:nTime
            dirSun = dirSunMat(iTime, :)';
            coverageMat(iIn, iTime) = computeCoverage(rObsMat, rTrgMat, dirSun, R, deltaR, CoverageOptions);
        end
        disp("Inclination: "+iIn+" / "+nIn+", time: "+toc)
    end
    coverageMat = coverageMat * 100;
    figure()
    plot(rad2deg(inVec), max(coverageMat,[],2) - min(coverageMat,[],2))
    xlabel("inclination [deg]")
    ylabel("coverage variation [%]")
    title("Triple coverage yearly variation vs inclination")
end

% figure()
% hold on
% plot3(dirSunMat(:,1), dirSunMat(:,2), dirSunMat(:,3))
% axis equal

%% Functions

% Simple function for propagating Sun direction (circular orbit assumption)
function dirSunMat = propagateDirSun(timeVec, dirSun0, Constants)
    % ECL   : ecliptic frame (inertial, ecliptic plane)
    % HCI   : heliocentric inertial frame (inertial, equatorial plane)
    % ECI   : Earth-centred inertial frame

    th = - Constants.ECL; % Ecliptic angle (Earth-Sun is -23.4°)
    N_t = size(timeVec, 2);

    % DCM rotating ECL to HCI
    dcmEclToHci = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];

    dirSunMat = nan(N_t, 3);
    dirEarthEcl0 = - dcmEclToHci' * dirSun0;
    if abs(dirEarthEcl0(3)) > 1e-5
        warning("Sun does not lie on ecliptic.")
    end
    theta0 = atan2(dirEarthEcl0(2),dirEarthEcl0(1));
    for k = 1:N_t
        % Earth direction in ECL
        dirEarthEclX = cos(theta0 + Constants.OM_S*timeVec(k));
        dirEarthEclY = sin(theta0 + Constants.OM_S*timeVec(k));
        dirEarthEcl = [dirEarthEclX dirEarthEclY 0]';

        % Sun direction in ECI
         dirSunTemp = - dcmEclToHci * dirEarthEcl;
         dirSunMat(k,:) = dirSunTemp';
    end
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