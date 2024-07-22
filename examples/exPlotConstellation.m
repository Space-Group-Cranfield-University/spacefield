% This script lets a user visualise a Walker constellation of given
% parameters. The default values are the result of the optimisation of a
% constellation of 400 observer satellites and 460 km altitude. See
% exConstellationOptimisation.m.

clc
clear
close all

%% Settings

Constants = initialiseAstronomicalConstants();
Parameters.nOrb = 8;
Parameters.nSatOrb = 50;
Parameters.h = 460;
Parameters.in = deg2rad(70);
Parameters.raan = 0.0625;
Parameters.dtheta = 0.0951;
sizeThetaVec = 100;
nSat = Parameters.nOrb * Parameters.nSatOrb;

%% Initialise constellation

Constellation = initialiseConstellation(Parameters, Constants);
Orbits = computeOrbits(sizeThetaVec, Parameters, Constants);

%% Collect all satellite points (otherwise plot is slow)

rSatMat = nan(nSat, 3);
for iSat = 1:nSat
    rSatMat(iSat,:) = Constellation(iSat).x0(1:3,1)';
end

%% Plot constellation

[X, Y, Z] = getEarthShape(Constants);

figure()
axis equal
set(gca,'Color','k')
hold on
xlabel("X, ECI [km]")
ylabel("Y, ECI [km]")
zlabel("Z, ECI [km]")

% Earth surface
surf(X, Y, Z, "EdgeColor", "none", "FaceColor", "b")

% Orbits
for iOrbits = 1:Parameters.nOrb
    plot3(Orbits(iOrbits).x(:,1), Orbits(iOrbits).x(:,2), Orbits(iOrbits).x(:,3), "w")
end

% Satellites
scatter3(rSatMat(:,1), rSatMat(:,2), rSatMat(:,3), "p", "filled")

%% Functions

function Orbits = computeOrbits(sizeThetaVec, Parameters, Constants)
    thetaVec = linspace(0, 2*pi, sizeThetaVec);
    for iOrbits = 1:Parameters.nOrb
        for iTheta = 1:sizeThetaVec
            kepTemp = [Constants.R_E + Parameters.h, 1e-6, Parameters.in, ...
                    Parameters.raan + (iOrbits - 1) * 2*pi / Parameters.nOrb, ...
                    0, thetaVec(iTheta)];
            Orbits(iOrbits).x(iTheta, :) = convertKepToCart(kepTemp, Constants.MU_E);
        end
    end
end

function [X, Y, Z] = getEarthShape(Constants)
    [X, Y, Z] = sphere;
    X = X*Constants.R_E;
    Y = Y*Constants.R_E;
    Z = Z*Constants.R_E;
end