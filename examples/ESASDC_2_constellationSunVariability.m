% This script analyses how coverage varies with three different
% constellations for optical surveillance of space debris. In particular,
% the script plots the variation of the average daily coverage over one
% year.

clear
clc
close all

%% Constants
Constants = initialiseAstronomicalConstants;
%Constants.ECL = 0;

%% Parameters
nTrg = 1e4; % Number of targets for Monte Carlo evaluation
MapParameters.sizeMap = 128; % Vertical size of coverage map (nTrgMap = 2 * sizeMap^2)
MapParameters.trgSizeMap = 0.7; % Size for generating coverage map. -1 sets to the average size above 10 cm (about 70 cm)
MapParameters.altitudeMap = 800; % Altitude for generating coverage map (800 km is the peak, 1100 is the average)
R = Constants.R_E + 100; % Horizon radius
dirSun0 = [-1 0 0]';

VisibilityParams.limitingMagnitude = +12;
VisibilityParams.F_spec     = @(x) 1/(4); % Specular scattering function
VisibilityParams.F_diff     = @(x) 2/(3*pi)*((pi-x)*cos(x)+sin(x)); % Diffusive scattering function
VisibilityParams.rho        = 0.2; % Average target reflectivity
VisibilityParams.beta       = 0.5; % beta = 1 scattering is completely diffusive, beta = 0 completely specular
VisibilityParams.Constants = Constants;

nOrbVec = [1 2 5];
nSatOrbVec = [20 10 4];
inVec = [97.3 97.3 75]; % 97.3 is SSO at 460 km
ConstellationParameters.h = 460;
ConstellationParameters.raan = deg2rad(90);
ConstellationParameters.phase_factor = 1; % up to N_sat
T = computeOrbitalPeriod(Constants.R_E + ConstellationParameters.h, Constants.MU_E);
n = sqrt(Constants.MU_E / (Constants.R_E + ConstellationParameters.h)^3);
deltaT = T / 4;
tspan = 0:deltaT:(T-deltaT);
MapParameters.timestepMap = 1;
tspanSun = linspace(0, 365*86400, 50);

legendConstellationType = ["single SSO","double SSO","Walker"];

%% Import data
importDirectory = input("Input data directory: ","s");
importAltitudeData
importDeclinationData
importSizeData

%% Preprocess

% Normalise data
AltitudeData{:,2} = AltitudeData{:,2}/sum(AltitudeData{:,2});
DeclinationData{:,2} = DeclinationData{:,2}/sum(DeclinationData{:,2});
SizeData{:,2} = SizeData{:,2}/sum(SizeData{:,2});

% Generate sample
sample = generateSample(nTrg, SizeData, AltitudeData, DeclinationData, Constants.R_E);
if MapParameters.trgSizeMap == -1
    MapParameters.trgSizeMap = mean(sample(:,4));
end

figure()
subplot(3,1,1)
plot(AltitudeData{:,1},AltitudeData{:,2}, "LineWidth", 2)
xlabel("altitude [km]", "FontSize", 12)
ylabel("density", "FontSize", 12)
fontsize(gca, 12, 'points')
subplot(3,1,2)
plot(DeclinationData{:,1}, DeclinationData{:,2}, "LineWidth", 2)
xlim([-90, 90])
xlabel("declination [deg]", "FontSize", 12)
ylabel("density", "FontSize", 12)
fontsize(gca, 12, 'points')
subplot(3,1,3)
semilogx(SizeData{:,1}, SizeData{:,2}, "LineWidth", 2)
xlabel("diameter [m]", "FontSize", 12)
ylabel("density", "FontSize", 12)
fontsize(gca, 12, 'points')

%% RAAN variation due to J2 effect
% Check https://control.asu.edu/Classes/MAE462/462Lecture13.pdf

raanDot = @(in) -1.5*n*Constants.J2*(Constants.R_E/(Constants.R_E + ConstellationParameters.h))^2*cos(in);
inVecTemp = deg2rad(0:0.1:180);
figure()
semilogy(rad2deg(inVecTemp), 2*pi./(abs(raanDot(inVecTemp)-Constants.OM_S)*86400))
fontsize(gca, 10, 'points')
ylabel("repeat period with respect to the Sun [days]", "LineWidth", 2)
xlabel("inclination [deg]", "LineWidth", 2)

%% Test Sun variability for multiple constellations

% Propagate Sun
dirSunMat = propagateDirSun(tspanSun, dirSun0, Constants);
coverageVec = [];
tripleCoverageVec = [];
for k = 1:size(legendConstellationType, 2)
    disp("")
    disp("Constellation: "+legendConstellationType(k))
    ConstellationParameters.in = deg2rad(inVec(k));
    ConstellationParameters.nOrb = nOrbVec(k);
    ConstellationParameters.nSatOrb = nSatOrbVec(k);
    ConstellationParameters.dtheta = ConstellationParameters.phase_factor*2*pi/(ConstellationParameters.nOrb * ConstellationParameters.nSatOrb);
    [coverageVecTemp, tripleCoverageVecTemp]     = mainCoverage(tspan, tspanSun,...
                raanDot, sample, dirSunMat, R, ...
                ConstellationParameters, VisibilityParams, Constants);
    coverageVec = [coverageVec coverageVecTemp'];
    tripleCoverageVec = [tripleCoverageVec tripleCoverageVecTemp'];
end

%% Plot

figure()
subplot(2,1,1)
hold on
for k = 1:3
    plot(tspanSun/86400, coverageVec(:, k)*100, "LineWidth", 2)
end
fontsize(gca, 10, 'points')
ylabel("single coverage [%]", "FontSize", 12)
xlabel("time [days]", "FontSize", 12)
xlim([0 365])

subplot(2,1,2)
hold on
for k = 1:3
    plot(tspanSun/86400, tripleCoverageVec(:, k)*100, "LineWidth", 2)
end
lgd = legend(legendConstellationType);
lgd.Location = "southwest";
fontsize(gca, 10, 'points')
ylabel("triple coverage [%]", "FontSize", 12)
xlabel("time [days]", "FontSize", 12)
xlim([0 365])

%% Function

function [coverageVec, tripleCoverageVec]     = mainCoverage(tspan, tspanSun,...
                raanDot, sample, dirSunMat, R, ...
                ConstellationParameters, VisibilityParameters, Constants)
   
    coverageVec = zeros(size(tspanSun));
    tripleCoverageVec = zeros(size(tspanSun));
    tic
    % 3) Compute coverage
    for k = 1:size(tspanSun, 2)
        ConstellationParameters.raan = ConstellationParameters.raan + raanDot(ConstellationParameters.in)*(tspanSun(2)-tspanSun(1));
        % 1) Initialise constellation
        Constellation = initialiseConstellation(ConstellationParameters, Constants);
        % 2) Propagate constellation
        Constellation = propagateConstellation(tspan, Constellation, Constants.MU_E);
        %disp("Propagation done, time: "+toc)
        coverageVec(k) = mean(computeCoverage(Constellation, 1, sample, dirSunMat(k,:)', R, VisibilityParameters));
        tripleCoverageVec(k) = mean(computeCoverage(Constellation, 3, sample, dirSunMat(k,:)', R, VisibilityParameters));
    end
    disp("Coverage computation done, time: "+toc)
end

function coverageVec = computeCoverage(Constellation, coverageFold, sample, dirSun, R, VisibilityParams)
    coverageVec = zeros(size(Constellation(1).x, 1), 1);
    for kT = 1:size(Constellation(1).x, 1)
        rObsMat = convertStructToMatAtTime(Constellation, kT);
        coverageVec(kT) = integrateMonteCarlo(sample, ...
                        @(sample) isTargetCovered(sample, rObsMat, dirSun, ...
                        R, VisibilityParams, coverageFold));
    end
end

function bool = isTargetCovered(sample, rObsMat, dirSun, R, VisibilityParams, coverageFold)
    count = 0;
    k = 1;
    while count < coverageFold && k < size(rObsMat, 1)
        rObs = rObsMat(k,:)';
        count = count + isVisible(rObs, sample(1:3), dirSun, R, sample(4), VisibilityParams);
        k = k+1;
    end
    bool = (count == coverageFold);
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

% Define visibility function
function flag = isVisible(rObs, rTrg, dirSun, R, trgSize, VisibilityParams)
    flag = isEclipsed(rTrg, dirSun, R);
    if flag == 0
        flag = ~isAboveTheHorizon(rObs, rTrg, R);
        if flag == 0
            flag = isScatteringWeakMag(rObs, rTrg, dirSun, trgSize, VisibilityParams);
        end
    end
    flag = ~flag;
end

% Counts how many sensors are available to the given target
function visibilityCount = countVisibility(rTrg, rObsMat, dirSun, R, trgSize, VisibilityParams)
    nObs = size(rObsMat, 1);
    visibilityCount = 0;
    iObs = 1;
    while iObs < nObs
        visibilityCount = visibilityCount + isVisible(rObsMat(iObs,:)', rTrg, dirSun, R, trgSize, VisibilityParams);
        iObs = iObs + 1;
    end
end

% Initialises a grid of target points on the surface of a sphere
function [rTrgMat, raDecTrgMat] = initialiseTargetPoints(nPoints, h, Constants)
    R = Constants.R_E + h;
    azVec = linspace(-pi, pi, 2*nPoints);
    elVec = linspace(pi/2, -pi/2, nPoints);
    rTrgMat = [];
    raDecTrgMat = [];
    for kAz = 1:(2*nPoints)
        az = azVec(kAz);
        for kEl = 1:(nPoints)
            el = elVec(kEl);
            x = R*cos(el)*cos(az);
            y = R*cos(el)*sin(az);
            z = R*sin(el);
            rTrg = [x y z];
            raDecTrg = getRaDec(rTrg);
            rTrgMat = [rTrgMat; rTrg];
            raDecTrgMat = [raDecTrgMat; raDecTrg];
        end
    end
end

% Counts visibility for a group of targets
function visibilityCountVec = countVisibilityMatrix(rTrgMat, rObsMat, dirSun, R, trgSize, VisibilityParams)
    nTrg = size(rTrgMat, 1);
    visibilityCountVec = zeros(1, nTrg);
    for iTrg = 1:nTrg
        visibilityCountVec(iTrg) = countVisibility(rTrgMat(iTrg, :)', rObsMat, dirSun, R, trgSize, VisibilityParams);
    end
end

% Converts vec to mat formatting for surf, with size(az)*size(el) =
% (2*nPoints)*(nPoints)
function [raMat, decMat, visibilityCountMat] = convertVecToMat(raDecTrgMat, visibilityCountVec, nPoints)
    raMat = raDecTrgMat(1:(nPoints), 1);
    decMat = raDecTrgMat(1:(nPoints), 2);
    visibilityCountMat = visibilityCountVec(1, 1:(nPoints))';
    for k = 2:(2*nPoints)
        iStart = (k-1)*(nPoints)+1;
        iEnd = k*(nPoints);
        raMat = [raMat raDecTrgMat(iStart:iEnd, 1)];
        decMat = [decMat raDecTrgMat(iStart:iEnd, 2)];
        visibilityCountMat = [visibilityCountMat visibilityCountVec(1, iStart:iEnd)'];
    end
    raMat = raMat(1:end-1, 1:end-1);
    decMat = decMat(1:end-1, 1:end-1);
    visibilityCountMat = visibilityCountMat(1:end-1, 1:end-1);
end

function magnitude = computeApparentMagnitude(rObs, rTrg, dirSun, trgSize, rho, F_diff, F_spec, beta, Constants)
    % NOTICE: this assumes target and observer are approximately at the
    % same distance from the Sun. Does not hold for arbitrary astronomical
    % object, only for targets that are debris in Earth orbit!
    % Inputs:
    %   - rObs      : vector containing observer position
    %   - rTrg      : vector containing target position
    %   - dirSun    : unit vector for Sun direction
    %   - trgSize   : size of target in m
    %   - rho       : target reflectivity
    %   - F_diff    : diffusive scattering as a function of phase angle
    %   - F_spec    : specular scattering as a function of phase angle
    %   - beta      : diffusive/specular mixing (can be between 0 and 1)
    %   - Constants : astronomical constants

    dr = rObs - rTrg;
    deltaR = norm(dr);
    phi = acos(dr' * dirSun/deltaR);
    % Compute magnitude
    V = Constants.M_S - 2.5*log10( rho * (trgSize*1e-3/2) ^ 2 / Constants.R_S ^ 2 ); % absolute magnitude target
    magnitude   = V + 5*log10( deltaR / Constants.R_S ) ...
                - 2.5 * log10( beta * F_diff(phi) + ( 1 - beta ) * F_spec(phi) ); % apparent magnitude target
end

function bool = isScatteringWeakMag(rObs, rTrg, dirSun, trgSize, VisibilityParams)
    % Scattering function based on limiting magnitude
    rho     = VisibilityParams.rho;
    F_diff  = VisibilityParams.F_diff;
    F_spec  = VisibilityParams.F_spec;
    beta    = VisibilityParams.beta;
    lMag    = VisibilityParams.limitingMagnitude;

    mag = computeApparentMagnitude(rObs, rTrg, dirSun, trgSize, rho, F_diff, F_spec, beta, VisibilityParams.Constants);
    bool = (mag > lMag);
end

% Transform constellation structure in matrix
function rObsMat = convertStructToMat(Constellation)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x0(1:3,1);
    end
end

% Transform constellation structure in matrix
function rObsMat = convertStructToMatAtTime(Constellation, timeStep)
    nConstellation = size(Constellation, 2);
    rObsMat = nan(nConstellation, 3);
    for iConstellation = 1:nConstellation
        rObsMat(iConstellation, :) = Constellation(iConstellation).x(timeStep, 1:3)';
    end
end

% gets right ascension and declination of a Cartesian vector
function raDec = getRaDec(r)
    ra = atan2(r(2), r(1));
    dec = asin(r(3)/norm(r));
    raDec = [ra dec];
end

% propagates constellation orbits and computes their radec
function Orbits = propagateOrbits(sizeThetaVec, Parameters, Constants)
    thetaVec = linspace(0, 2*pi, sizeThetaVec);
    for iOrbits = 1:Parameters.nOrb
        for iTheta = 1:sizeThetaVec
            kepTemp = [Constants.R_E + Parameters.h, 1e-6, Parameters.in, ...
                    Parameters.raan + (iOrbits - 1) * 2*pi / Parameters.nOrb, ...
                    0, thetaVec(iTheta)];
            Orbits(iOrbits).x(iTheta, :) = convertKepToCart(kepTemp, Constants.MU_E);
            Orbits(iOrbits).raDec(iTheta, :) = getRaDec(Orbits(iOrbits).x(iTheta, :));
        end
    end
end

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

function sample = generateSample(nTrg, SizeData, AltitudeData, DeclinationData, R_E)
    % Function for generating samples to feed to the Monte Carlo
    % integrator. A sample contains the three Cartesian coordinates and the
    % debris size specification. The number of samples is nTrg, indicating
    % the number of virtual targets obeserved by the system.

    sizeSample = randsample(SizeData{:,1}, nTrg, true, SizeData{:,2});
    altitudeSample = randsample(AltitudeData{:,1}, nTrg, true, AltitudeData{:,2});
    declinationSample = deg2rad(randsample(DeclinationData{:,1}, nTrg, true, DeclinationData{:,2}));
    sample = zeros(nTrg, 4);
    for k = 1:nTrg
        raSampleTemp = rand()*2*pi;
        sizeSampleTemp = sizeSample(k);
        altitudeSampleTemp = altitudeSample(k);
        declinationSampleTemp = declinationSample(k);
        radiusTemp = altitudeSampleTemp + R_E;
        x = radiusTemp * cos(raSampleTemp) * cos(declinationSampleTemp);
        y = radiusTemp * sin(raSampleTemp) * cos(declinationSampleTemp);
        z = radiusTemp * sin(declinationSampleTemp);
        sampleTemp = [x, y, z, sizeSampleTemp];
        sample(k, :) = sampleTemp;
    end
end