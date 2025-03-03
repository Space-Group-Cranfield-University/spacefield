% Simulates the distribution of relative velocities for random targets in
% circular orbits (LEO).
%
% This script replicates results shown in Fig. 5 at https://doi.org/10.1016/j.actaastro.2025.02.019

clear
close all
clc

% Inputs
Constants = initialiseAstronomicalConstants;
altitudeObserver = 400; % [km]
phi = deg2rad(90); % target inclination
deltaR = 10; % relative distance [km]
N_MC = 1e5;
N_PHI = 50;
N_R = 50;
phiVec = linspace(0,deg2rad(90),N_PHI);
deltaRVec = 10.^linspace(0, 3.5, N_R); 
positionNormObserver = Constants.R_E + altitudeObserver;
velocityNormObserver = computeVelocityNorm(positionNormObserver, Constants);
epsilon = deltaR/positionNormObserver;
test = "tangential velocity vs r";

% 1) Compute relative velocity distribution at phi = 0

dvVec_1 = zeros(1,N_MC);
phi = deg2rad(90);
for k = 1:N_MC
    dvVec_1(k) = computeRelativeVelocityFixedOrientation(altitudeObserver, deltaR, phi, Constants);
end

figure()
histogram(dvVec_1)

% 2) Compute expected value relative velocity vs phi

if strcmp(test, "relative velocity vs phi")
    dvMat_2 = zeros(N_PHI, N_MC);
    for kPhi = 1:N_PHI
        for kMC = 1:N_MC
            dvMat_2(kPhi, kMC) = computeRelativeVelocityFixedOrientation(altitudeObserver, deltaR, phiVec(kPhi), Constants);
        end
        disp(kPhi)
    end
    dvAvg_2 = mean(dvMat_2, 2);
    expectedValueFunction_2 = @(phi) 2*sqrt(2)/pi * velocityNormObserver * sqrt(2+cos(phi)*epsilon);
    
    figure()
    plot(rad2deg(phiVec), dvAvg_2)
    hold on
    plot(rad2deg(phiVec), expectedValueFunction_2(phiVec))
end

% 3) Compute expected value relative velocity vs deltaR

if strcmp(test, "relative velocity vs r")
    dvMat_3 = zeros(N_R, N_MC);
    expectedValue_3 = zeros(1, N_R);
    for kR = 1:N_R
        for kMC = 1:N_MC
            dvMat_3(kR, kMC) = computeRelativeVelocityVariableOrientation(altitudeObserver, deltaRVec(kR), Constants);
        end
        integrand = @(phi) 2/pi*sqrt(2 + cos(phi)*deltaRVec(kR)/positionNormObserver);
        compIntegral = integral(integrand, 0, pi/2);
        expectedValue_3(kR) = compIntegral * 2 * sqrt(2) / pi * velocityNormObserver;
        disp(kR)
    end
    
    dvAvg_3 = mean(dvMat_3,2);
    
    figure()
    semilogx(deltaRVec, dvAvg_3)
    hold on
    semilogx(deltaRVec, expectedValue_3)
end

% 4) Compute expected value tangential velocity vs phi

if strcmp(test, "tangential velocity vs phi")    
    vTMat_4 = zeros(N_PHI, N_MC);
    for kPhi = 1:N_PHI
        for kMC = 1:N_MC
            vTMat_4(kPhi, kMC) = computeTangentialVelocityFixedOrientation(altitudeObserver, deltaR, phiVec(kPhi), Constants);
        end
        disp(kPhi)
    end
    vTAvg_4 = mean(vTMat_4, 2);
    analyticalResults(1) = sqrt(2) * 2 / pi * sqrt(2) * velocityNormObserver;
    analyticalResults(2) = 8 / pi^2 * velocityNormObserver;

    figure()
    plot(rad2deg(phiVec), vTAvg_4);
    hold on
    plot([rad2deg(phiVec(1)) rad2deg(phiVec(end))], analyticalResults)
end

% 5) Compute expected value tangential velocity vs dr

if strcmp(test, "tangential velocity vs r")    
    vTMat_5 = zeros(N_R, N_MC);
    deltaVMat_5 = zeros(N_R, N_MC);
    for kR = 1:N_R
        for kMC = 1:N_MC
            vTMat_5(kR, kMC) = computeTangentialVelocityVariableOrientation(altitudeObserver, deltaRVec(kR), Constants);
            %deltaVMat_5(kR, kMC) = computeRelativeVelocityVariableOrientation(altitudeObserver, deltaRVec(kR), Constants);
        end
        disp(kR)
    end
    vTAvg_5 = mean(vTMat_5, 2);
    %deltaVAvg_5 = mean(deltaVMat_5, 2);
    analyticalResults = @(r) velocityNormObserver * sqrt(2 + r/positionNormObserver) * sqrt(2) * 2 / pi;

    figure()

    subplot(2,1,1)
    semilogx(deltaRVec, vTAvg_5, "LineWidth", 2);
    hold on
    %semilogx(deltaRVec, deltaVAvg_5, "LineWidth", 2)
    semilogx(deltaRVec, analyticalResults(deltaRVec), "LineWidth", 2)
    fontsize(gca, 12, 'points')
    legend(["Monte Carlo", "Analytical"])
    %title("Expected tangential velocity", "FontSize", 12)
    xlabel("relative distance [km]", "FontSize", 12)
    ylabel("expected velocity [km/s]", "FontSize", 12)

    subplot(2,1,2)
    semilogx(deltaRVec, abs(100*(vTAvg_5' - analyticalResults(deltaRVec)))./vTAvg_5', "LineWidth", 2)
    fontsize(gca, 12, 'points')
    %title("Numerical-Analytical relative error", "FontSize", 12)
    xlabel("relative distance [km]", "FontSize", 12)
    ylabel("relative error [%]", "FontSize", 12)
end

%% Functions

function v = computeVelocityNorm(a, Constants)
    v = sqrt(Constants.MU_E / a);    
end

function deltaV = computeRelativeVelocityFixedOrientation(altitudeObserver, deltaR, phi, Constants)
    % Compute observer position
    positionNormObserver = Constants.R_E + altitudeObserver;
    velocityNormObserver = computeVelocityNorm(positionNormObserver, Constants);
    positionObserver = [positionNormObserver 0 0]';
    
    % Random observer velocity
    psi = rand() * 2*pi;
    velocityObserver = velocityNormObserver * [0 cos(psi) sin(psi)]';
    
    % Compute target position
    deltaPositionTarget = deltaR * [cos(phi) sin(phi) 0]';
    positionTarget = positionObserver + deltaPositionTarget;
    positionNormTarget = norm(positionTarget);
    
    % Random target velocity
    velocityNormTarget = computeVelocityNorm(positionNormTarget, Constants);
    theta = rand() * 2*pi;
    velocityTargetBeforeRotation = velocityNormTarget * [0 cos(theta) sin(theta)]';
    %alpha = atan(sin(phi) / (cos(phi) + positionNormObserver / deltaR));
    RotationMatrix = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    velocityTarget = RotationMatrix * velocityTargetBeforeRotation;
    
    % Compute relative velocity
    n = cross(positionObserver, velocityObserver)/positionNormObserver^2;
    deltaVelocityVector = velocityTarget - velocityObserver - cross(n, deltaPositionTarget);
    deltaV = norm(deltaVelocityVector);
end

function deltaV = computeRelativeVelocityVariableOrientation(altitudeObserver, deltaR, Constants)
    phi = rand() * pi/2;
    deltaV = computeRelativeVelocityFixedOrientation(altitudeObserver, deltaR, phi, Constants);
end



function vT = computeTangentialVelocityFixedOrientation(altitudeObserver, deltaR, phi, Constants)
    % Compute observer position
    positionNormObserver = Constants.R_E + altitudeObserver;
    velocityNormObserver = computeVelocityNorm(positionNormObserver, Constants);
    positionObserver = [positionNormObserver 0 0]';
    
    % Random observer velocity
    psi = rand() * 2*pi;
    velocityObserver = velocityNormObserver * [0 cos(psi) sin(psi)]';
    
    % Compute target position
    deltaPositionDirTarget = [cos(phi) sin(phi) 0]';
    deltaPositionTarget = deltaR * deltaPositionDirTarget;
    positionTarget = positionObserver + deltaPositionTarget;
    positionNormTarget = norm(positionTarget);
    
    % Random target velocity
    velocityNormTarget = computeVelocityNorm(positionNormTarget, Constants);
    theta = rand() * 2*pi;
    velocityTargetBeforeRotation = velocityNormTarget * [0 cos(theta) sin(theta)]';
    %alpha = atan(sin(phi) / (cos(phi) + positionNormObserver / deltaR));
    RotationMatrix = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    velocityTarget = RotationMatrix * velocityTargetBeforeRotation;
    
    % Compute relative velocity
    n = cross(positionObserver, velocityObserver)/positionNormObserver^2;
    deltaVelocityVector = velocityTarget - velocityObserver - cross(n, deltaPositionTarget);
    %deltaV = norm(deltaVelocityVector);

    % Compute tangential velocity
    vTVec = deltaVelocityVector - dot(deltaVelocityVector, deltaPositionDirTarget)*deltaPositionDirTarget;
    vT = norm(vTVec);
end

function vT = computeTangentialVelocityVariableOrientation(altitudeObserver, deltaR, Constants)
    phi = rand() * pi/2;
    vT = computeTangentialVelocityFixedOrientation(altitudeObserver, deltaR, phi, Constants);
end