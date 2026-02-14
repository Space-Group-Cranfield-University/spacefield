% This script analyses how the yearly Sun position affects the coverage of
% a given constellation, both for fixed and variable constellation orbit
% inclination.

% This script replicates results shown in Figures 15-16 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
clc
close all

Constants = initialiseAstronomicalConstants;
dirSun0 = [-1 0 0]';
%dirSun1 = [0 -cosd(Constants.ECL) sind(Constants.ECL)]';
nTime = 20; % substitute 40
timeVec = linspace(0, 52*7*86400, nTime);
ConstellationParameters.nOrb = 5;
ConstellationParameters.nSatOrb = 4;
ConstellationParameters.h = 460;
ConstellationParameters.in = deg2rad(70);
ConstellationParameters.raan = 0;
ConstellationParameters.dtheta = deg2rad(5);
CoverageOptions.coverageCount = 1; % Triple coverage
R = Constants.R_E + 100;
deltaR = 2000;
SampleOptions.R_1 = Constants.R_E + 400;
SampleOptions.R_2 = Constants.R_E + 1800;
SampleOptions.N = 1e4; % substitute 5*1e4
inVec = deg2rad(linspace(65,85,9)); % vector of inclinations
nIn = size(inVec,2);
yearlyVariabilityAnalysis = 1;
parametricAnalysis = 0;

%% Sun propagation

dirSunMat = propagateDirSun(timeVec, dirSun0, Constants);
dirSunMatEclOnly = propagateDirSunEclOnly(timeVec, dirSun0, Constants);

%% Compute coverage

% Initialise Monte Carlo sample
rTrgMat = initialiseSampleShell(SampleOptions);

% Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);
rObsMat = convertStructToMat(Constellation);

% Coverage computation
if yearlyVariabilityAnalysis
    coverageVec = nan(1, nTime);
    coverageVecEclOnly = nan(1, nTime);
    tic
    for iTime = 1:nTime
        dirSun = dirSunMat(iTime, :)';
        dirSunEclOnly = dirSunMatEclOnly(iTime, :)';
        coverageVec(iTime)  = computeCoverage(rObsMat, rTrgMat, dirSun, ...
                            R, deltaR, CoverageOptions);
        coverageVecEclOnly(iTime)  = computeCoverage(rObsMat, rTrgMat, dirSunEclOnly, ...
                            R, deltaR, CoverageOptions);
        disp("Timestep: "+iTime+" / "+nTime+", time: "+toc)
    end  
    figure()

    subplot(2,1,1)
    plot(timeVec/86400, coverageVec*100, "LineWidth", 2)
    fontsize(gca, 12, 'points')
    %title("Yearly triple coverage", "FontSize", 12)
    xlabel("time [days]", "FontSize", 14)
    ylabel("coverage [%]", "FontSize", 14)
    xlim([0, 365])

    subplot(2,1,2)
    plot(timeVec/86400, coverageVecEclOnly*100, "LineWidth", 2)
    fontsize(gca, 12, 'points')
    %title("Yearly triple coverage, Sun elevation only", "FontSize", 12)
    xlabel("time [days]", "FontSize", 14)
    ylabel("coverage [%]", "FontSize", 14)
    xlim([0, 365])
    ylim([34, 37])

    figure()
    plot(timeVec/86400, coverageVec*100, "LineWidth", 2)
    hold on
    plot(timeVec/86400, coverageVecEclOnly*100, "LineWidth", 2)
    fontsize(gca, 12, 'points')
    title("Yearly triple coverage", "FontSize", 12)
    fontsize(gca, 12, 'points')
    xlabel("time [days]", "FontSize", 14)
    ylabel("coverage [%]", "FontSize", 14)
    xlim([0, 365])
    %ylim([34, 37])
    legend(["az-el","el only"])
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
            coverageMat(iIn, iTime) = computeCoverage(rObsMat, rTrgMat,...
                                    dirSun, R, deltaR, CoverageOptions);
        end
        disp("Inclination: "+iIn+" / "+nIn+", time: "+toc)
    end
    coverageMat = coverageMat * 100;

    figure()

    subplot(2,1,1)
    plot(rad2deg(inVec), mean(coverageMat,2), "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("inclination [deg]", "FontSize", 14)
    ylabel("avg. coverage [%]", "FontSize", 14)
    %title("Triple coverage: yearly average against Inclination", "FontSize", 12)
    xlim([65, 85])

    subplot(2,1,2)
    plot(rad2deg(inVec), max(coverageMat,[],2) - min(coverageMat,[],2), "LineWidth", 2)
    fontsize(gca, 12, 'points')
    xlabel("inclination [deg]", "FontSize", 14)
    ylabel("coverage var. [%]", "FontSize", 14)
    %title("Triple coverage: yearly variation against Inclination", "FontSize", 12)
    xlim([65, 85])
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

function dirSunMat = propagateDirSunEclOnly(timeVec, dirSun0, Constants)
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
        dirSunTemp = dcmEclToHci * dirEarthEcl;
        dirSunTemp(1:2,1) = dirSun0(1:2,1);
        dirSunTemp = dirSunTemp / norm(dirSunTemp);
        dirSunMat(k,:) = dirSunTemp';
    end
end

function dirSunMat = propagateDirSunAzOnly(timeVec, dirSun0, Constants)
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
        dirSunTemp = dcmEclToHci * dirEarthEcl;
        dirSunTemp(3,1) = dirSun0(3,1);
        dirSunTemp = dirSunTemp / norm(dirSunTemp);
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