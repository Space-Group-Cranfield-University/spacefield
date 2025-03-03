% This script lets a user visualise a Walker constellation of given
% parameters. The default values are the result of the optimisation of a
% constellation of 400 observer satellites and 460 km altitude. See
% exConstellationOptimisation.m.

% This script replicates results shown in Fig. 21 at 
% https://doi.org/10.1016/j.actaastro.2025.02.019 when the correct
% constellation parameters are substituted

clc
clear
close all

%% Settings

Constants = initialiseAstronomicalConstants();
Parameters.nOrb = 10;
Parameters.nSatOrb = 50;
Parameters.h = 500;
Parameters.in = deg2rad(70);
Parameters.raan = 0;
Parameters.nSatOrbActual = 5;
areOrbitsOverlapped = 0;
%Parameters.dtheta = deg2rad(8.5) - areOrbitsOverlapped * 2*pi * Parameters.nSatOrbActual / Parameters.nSatOrb;
Parameters.dtheta = deg2rad(10);
sizeThetaVec = 100;
dirSun = [-1 0 0]';
deltaR = 2000;
R = Constants.R_E + 100;

% Constellation projection
h_1 = 800; % Peak altitude 1
h_2 = 800; % Peak altitude 2
nPoints = 128;

%% Initialise constellation

%Constellation = initialiseLowBudgetConstellation(Parameters, Constants);
Constellation = initialiseConstellation(Parameters, Constants);
nSat = size(Constellation, 2);
Orbits = propagateOrbits(sizeThetaVec, Parameters, Constants);
rObsMat = convertStructToMat(Constellation);

%% Initialise target

Target.kep0Est  = [Constants.R_E + 800, 1e-5, deg2rad(98), deg2rad(90), 0, 0]';
tend = computeOrbitalPeriod(Target.kep0Est(1), Constants.MU_E);
tspan = linspace(0, tend, 100);
Target.x0Est = convertKepToCart(Target.kep0Est, Constants.MU_E); % Estimated initial state
Target.x = propagateKeplerian(tspan, Target.x0Est, Constants.MU_E);

%% Initialise target points

[rTrgMat_1, raDecTrgMat_1] = initialiseTargetPoints(nPoints, h_1, Constants);
[rTrgMat_2, raDecTrgMat_2] = initialiseTargetPoints(nPoints, h_2, Constants);

%% Compute coverage

visibilityCountVec_1 = countVisibilityMatrix(rTrgMat_1, rObsMat, dirSun, R, deltaR);
visibilityCountVec_2 = countVisibilityMatrix(rTrgMat_2, rObsMat, dirSun, R, deltaR);
[raMat_1, decMat_1, visibilityCountMat_1] = convertVecToMat(raDecTrgMat_1, visibilityCountVec_1, nPoints);
[raMat_2, decMat_2, visibilityCountMat_2] = convertVecToMat(raDecTrgMat_2, visibilityCountVec_2, nPoints);

%% Collect all satellite points (otherwise plot is slow)

rSatMat = nan(nSat, 3);
for iSat = 1:nSat
    rSatMat(iSat,:) = Constellation(iSat).x0(1:3,1)';
end

%% Plot constellation

% 3D plot options
Az = 135;
El = 30;
figure()
plotConstellation(Parameters, Orbits, rSatMat, Az, El, 1, 1)
hold on
plot3(Target.x(:,1), Target.x(:,2), Target.x(:,3), "c", "LineWidth", 2)

% 3D plot options
Az = 135;
El = 30;
figure()
plotConstellation(Parameters, Orbits, rSatMat, Az, El, 1, 1)

%% Plot radec

% figure()
% hold on
% for iOrbits = 1:Parameters.nOrb
%     plotRaDec(Orbits(iOrbits).raDec, "#EDB120", 1)
% end
% axis equal
% ylim([-90, +90])
% xlim([-180, 180])

%% Plot visibility map

figure()

subplot(2,1,1)
m = mesh(rad2deg(raMat_1), rad2deg(decMat_1), (visibilityCountMat_1>2));
m.FaceColor = 'flat';
view(2)
axis equal
cb = colorbar;
ylim([-90, +90])
xlim([-180, 180])
hold on
for iOrbits = 1:Parameters.nOrb
    plotRaDec3D(Orbits(iOrbits).raDec, "red", 1, 1e5)
end
ylabel("declination [deg]",'FontSize',14)
ylabel(cb,'triple visibility','FontSize',14,'Rotation',270)

subplot(2,1,2)
m = mesh(rad2deg(raMat_1), rad2deg(decMat_1), visibilityCountMat_1);
m.FaceColor = 'flat';
view(2)
axis equal
cb = colorbar;
%clim([0, 8])
ylim([-90, +90])
xlim([-180, 180])
ylabel("declination [deg]",'FontSize',14)
xlabel("right ascension [deg]",'FontSize',14)
ylabel(cb,'visibility count','FontSize',14,'Rotation',270)

%% Functions

% Keplerian propagation function (without ode)
function xMat = propagateKeplerian(tspan, x0, mu)
    kep0 = convertCart2Kep(x0, mu);
    n = sqrt(mu / kep0(1)^3);
    xMat = nan(size(tspan, 2), 6);
    for k = 1:size(tspan, 2)
        dt = tspan(k) - tspan(1);
        dM = n * dt;
        dv = convertAnomalyM2v(dM, kep0(2));
        kepCurrent = [kep0(1:5); kep0(6) + dv];
        xMat(k, :) =  convertKepToCart(kepCurrent, mu)';
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
    flag    = ~( (~isWithinRange(rObs, rTrg, deltaR)) || ... 
            (isScatteringWeak(rObs, rTrg, dirSun)) || (isEclipsed(rTrg, dirSun, R)) || ...
            (~isAboveTheHorizon(rObs, rTrg, R)) );
    %flag    = ~ isEclipsed(rTrg, dirSun, R);
end

% Counts visibility for a group of targets
function visibilityCountVec = countVisibilityMatrix(rTrgMat, rObsMat, dirSun, R, deltaR)
    nTrg = size(rTrgMat, 1);
    visibilityCountVec = zeros(1, nTrg);
    for iTrg = 1:nTrg
        visibilityCountVec(iTrg) = countVisibility(rTrgMat(iTrg, :)', rObsMat, dirSun, R, deltaR);
    end
end

% Counts how many sensors are available to the given target
function visibilityCount = countVisibility(rTrg, rObsMat, dirSun, R, deltaR)
    nObs = size(rObsMat, 1);
    visibilityCount = 0;
    iObs = 1;
    while iObs < nObs
        visibilityCount = visibilityCount + isVisible(rObsMat(iObs,:)', rTrg, dirSun, R, deltaR);
        iObs = iObs + 1;
    end
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

% gets right ascension and declination of a Cartesian vector
function raDec = getRaDec(r)
    ra = atan2(r(2), r(1));
    dec = asin(r(3)/norm(r));
    raDec = [ra dec];
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

% plots constellation in 3D around the Earth
function plotConstellation(Parameters, Orbits, rSatMat, az, el, orbitFlag, satFlag)

    axis equal
    %set(gca,'Color','k')
    hold on
    xlabel("X, ECI [km]")
    ylabel("Y, ECI [km]")
    zlabel("Z, ECI [km]")
    
    % Earth surface
    cdata = imread("2k_earth_daymap_clouds.jpg");
    [x,y,z] = sphere(100);
    %z = 0.001*z;
    Re = 6371;
    earth = surf(Re*x, Re*y, -Re*z);
    set(earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none'); % "texturemap"
    %surf(X, Y, Z, "EdgeColor", "none", "FaceColor", "b")
    
    % Orbits
    if orbitFlag
        for iOrbits = 1:Parameters.nOrb
            plot3(Orbits(iOrbits).x(:,1), Orbits(iOrbits).x(:,2), Orbits(iOrbits).x(:,3), "Color", "#EDB120", "LineWidth", 1)
        end
    end
    
    % Satellites
    if satFlag
        scatter3(rSatMat(:,1), rSatMat(:,2), rSatMat(:,3), "filled", "v", "MarkerFaceColor", "#D95319")
    end
    
    %lighting phong
    view(az, el)
    %set(gca, 'Projection','perspective')
end

% plots ground track without joining right ascension from 2pi to 0
function [] = plotRaDec(raDecMat, color, lineWidth)
    hold on
    nRaDec = size(raDecMat,1);
    raDecPlot = raDecMat(1,:);
    for k = 2:nRaDec
        if abs(raDecMat(k,1) - raDecMat(k-1,1)) > pi
            plot(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), "Color", color, "LineWidth", lineWidth)
            raDecPlot = raDecMat(k,:);
        else
            raDecPlot = [raDecPlot; raDecMat(k,:)];
        end
    end
    plot(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), "Color", color, "LineWidth", lineWidth)
end

% plots ground track without joining right ascension from 2pi to 0 in 3D
function [] = plotRaDec3D(raDecMat, color, lineWidth, Z)
    hold on
    nRaDec = size(raDecMat,1);
    raDecPlot = raDecMat(1,:);
    for k = 2:nRaDec
        if abs(raDecMat(k,1) - raDecMat(k-1,1)) > pi
            plot3(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), Z*ones(size(raDecPlot, 1),1), "Color", color, "LineWidth", lineWidth)
            raDecPlot = raDecMat(k,:);
        else
            raDecPlot = [raDecPlot; raDecMat(k,:)];
        end
    end
    plot3(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), Z*ones(size(raDecPlot, 1),1), "Color", color, "LineWidth", lineWidth)
end

%% Ground track + Target track
