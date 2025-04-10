% This script analyses the quality of tracks for a given surveillance
% constellation, assuming a tracking strategy and a limited FOV. In
% particular, the script plots the percentage of targets that have had at
% least one track within half a orbital period of the given timestep for
% various tracking strategies.

% NOTICE: 
% - you could reorder targets based on maximum distance from
% constellation orbit (closest first, farthest last). This should better
% show how measurement opportunities increase for closer targets (an
% intuitive fact).
% - to avoid data dependence, all constellations should be tested on the
% same population data. Moreover, tests should be repeated for different
% target populations.
%
% TODO: Orbital parameter distribution of detected vs undetected objects

clear
close all
clc

%% Constants
Constants = initialiseAstronomicalConstants;

%% Parameters
nTrg = 500; % number of targets
deltaT = 30; % [s], time between consecutive measurements
nOrbitalPeriods = 5;
nBins = 20; % for plotting track results
loadTargets = 1;
saveTargets = 0;
importDirectory = input("Input data directory: ","s");

trackingStrategy = ["Radial", "Along-track", "Cross-track", "Greedy"];
%trackingStrategy = ["Radial", "Along-track", "Cross-track"];
%trackingStrategy = ["Greedy"];
ConstellationParameters.h = 460;
ConstellationParameters.raan = deg2rad(90);
ConstellationParameters.phase_factor = 1; % up to N_sat
ConstellationParameters.in = deg2rad(97.3);
ConstellationParameters.nOrb = 2;
ConstellationParameters.nSatOrb = 10;
dirSun = [-1 0 0]';
R = Constants.R_E + 100; % Horizon radius
VisibilityParams.limitingMagnitude = +12;
VisibilityParams.F_spec     = @(x) 1/(4); % Specular scattering function
%VisibilityParams.F_diff     = @(x) 2/(3*pi)*(x > pi/2); % Diffusive scattering function
VisibilityParams.F_diff     = @(x) 2/(3*pi)*((pi-x)*cos(x)+sin(x)); % Diffusive scattering function
VisibilityParams.rho        = 0.3; % Average target reflectivity
VisibilityParams.beta       = 1; % beta = 1 scattering is completely diffusive, beta = 0 completely specular
VisibilityParams.Constants = Constants;
VisibilityParams.FOV = deg2rad(5); % Sensor field of view
VisibilityParams.halfFOV = VisibilityParams.FOV / 2;
maxSlewRate = deg2rad(1); % max 1 degree/s
maxSlewRate = maxSlewRate * deltaT;
nRunsGreedy = 20;

nObs = ConstellationParameters.nOrb * ConstellationParameters.nSatOrb;
ConstellationParameters.dtheta = ConstellationParameters.phase_factor*2*pi/(ConstellationParameters.nOrb * ConstellationParameters.nSatOrb);
ConstellationParameters.T = computeOrbitalPeriod(Constants.R_E + ConstellationParameters.h, Constants.MU_E);
tspan = 1:deltaT:(nOrbitalPeriods*ConstellationParameters.T);
halfOrbitalPeriodLength = floor(ConstellationParameters.T/(2*deltaT));
timestepsBeforeTrackLoss = halfOrbitalPeriodLength;
minTrackLength = 10; % timesteps
minTrackSeparation = 10; % timesteps

%% Import and process population data
tic
if loadTargets
    %% 0-2) Import population
    load(importDirectory+"\Target.mat")
    disp("Targets loaded succesfully")
else
    %% 0) Import population data
    importLEOData
    importSizeData
    SizeData{:,2} = SizeData{:,2}/sum(SizeData{:,2});
    periodVec = zeros(size(SATCAT, 1),1);
    inclinationVec = zeros(size(SATCAT, 1),1);
    apogeeVec = zeros(size(SATCAT, 1),1);
    perigeeVec = zeros(size(SATCAT, 1),1);
    eccentricityVec = zeros(size(SATCAT, 1),1);
    smaVec = zeros(size(SATCAT, 1),1);
    for k = 1:size(SATCAT, 1)
        periodVec(k) = SATCAT.PERIOD(k);
        inclinationVec(k) = SATCAT.INCLINATION(k);
        apogeeVec(k) = SATCAT.APOGEE(k);
        perigeeVec(k) = SATCAT.PERIGEE(k);
        smaVec(k) = (apogeeVec(k) + perigeeVec(k)) / 2;
        eccentricityVec(k) = (apogeeVec(k) - perigeeVec(k)) / (apogeeVec(k) + perigeeVec(k));
    end
    smaVec = smaVec + Constants.R_E;
    disp("Data import completed, time: "+toc)
    
    %% 1) Generate sample target population
    for k = 1:nTrg
        eccTemp = 1;
        while eccTemp > 1e-1
            eccTemp = randsample(eccentricityVec, 1);
        end
        Target(k).kep   = [randsample(smaVec, 1), eccTemp, ...
                        deg2rad(randsample(inclinationVec, 1)), rand*2*pi, ...
                        rand*2*pi, rand*2*pi];
        Target(k).x   = convertKepToCart(Target(k).kep', Constants.MU_E)';
        Target(k).size  = randsample(SizeData{:,1}, 1, true, SizeData{:,2});
        Target(k).n     = sqrt(Constants.MU_E / Target(k).kep(1)^3);
        Target(k).p     = Target(k).kep(1) * (1-Target(k).kep(2)^2);
        Target(k).raanDot   = -1.5*Constants.J2*cos(Target(k).kep(3))*Target(k).n*(Constants.R_E / Target(k).p)^2;
        Target(k).aopDot    = 1.5*Constants.J2*(2-2.5*sin(Target(k).kep(3))^2)*Target(k).n*(Constants.R_E / Target(k).p)^2;
    end
    %% 2) Propagate targets (takes 3 mins for 500 targets, 4 s for 10 targets)
    for j = 2:size(tspan, 2)
        for k = 1:nTrg
            Target(k) = propagateJ2(deltaT, Target(k));
            Target(k).x(j, :) = convertKepToCart(Target(k).kep(j, :), Constants.MU_E);
        end
    end
    disp("Target propagation completed, time: "+toc)
    if saveTargets
        save(importDirectory+"\Target.mat", "Target")
    end
end
disp("Target initialisation completed, time: "+toc)

%% MAIN
nObs = ConstellationParameters.nOrb * ConstellationParameters.nSatOrb;
ConstellationParameters.dtheta = ConstellationParameters.phase_factor*2*pi/(ConstellationParameters.nOrb * ConstellationParameters.nSatOrb);

%% 3) Initialise constellation
Constellation = initialiseConstellation(ConstellationParameters, Constants);
disp("Constellation initialised, time: "+toc)

%% 4) Propagate constellation
Constellation = propagateConstellation(tspan, Constellation, Constants.MU_E);
disp("Constellation propagated, time: "+toc)

%% 4.5) Compute DCM_ECI2LVLH for constellation satellites at all times
for k = 1:size(Constellation, 2)
    for j = 1:size(tspan, 2)
        Constellation(k).DCM_ECI2LVLH(:,:,j) = get_DCM_ECI2LVLH(Constellation(k).x(j, :)');
    end
end
disp("DCM computed, time: "+toc)
disp(" ")

%% 5) Compute measurements (takes 30 s for 500 targets, 20 observers)
% This section is particuarly slow, mostly due to
% isTargetVisibleToObserver (both self and eclipsed, aboveTheHorizon)

numberOfPotentialMeasurements = zeros(size(tspan, 2), size(trackingStrategy,2));
nObservedTargetsHalfOrbitalPeriod = zeros(size(trackingStrategy, 2), size(tspan, 2));
nObservedTargets = zeros(size(trackingStrategy, 2), size(tspan, 2));
for kStrategy = 1:size(trackingStrategy, 2)
    disp(" ")
    disp("Tracking strategy: "+trackingStrategy(kStrategy))
    visibilityMat = zeros(nTrg, size(tspan, 2));
    for j = 1:size(tspan, 2)
        for kObs = 1:nObs
            rObs = Constellation(kObs).x(j, 1:3)';
            if strcmp(trackingStrategy(kStrategy), "Greedy") && j > 1
                Constellation(kObs).orientation(j, :)   = getOrientation(trackingStrategy(kStrategy), ...
                                                        Constellation(kObs).orientation(j-1, :)',...
                                                        j, Target, rObs, ...
                                                        Constellation(kObs).DCM_ECI2LVLH(:, :, j), ...
                                                        dirSun, R, Target(kTrg).size, VisibilityParams, ...
                                                        maxSlewRate, nRunsGreedy)';
            else
                Constellation(kObs).orientation(j, :) = getOrientation(trackingStrategy(kStrategy))';
            end
            Constellation(kObs).dirOrientation(j, :) = getDirOrientation(Constellation(kObs).orientation(j, :))';
            for kTrg = 1:nTrg
                rTrg = Target(kTrg).x(j, 1:3)';
                flag    = isTargetVisibleToObserver(rTrg, rObs, ...
                        Constellation(kObs).dirOrientation(j, :)', ...
                        Constellation(kObs).DCM_ECI2LVLH(:, :, j), dirSun, R, ...
                        Target(kTrg).size, VisibilityParams);
                numberOfPotentialMeasurements(j, kStrategy) = numberOfPotentialMeasurements(j, kStrategy) + flag;
                visibilityMat(kTrg, j) = visibilityMat(kTrg, j) + flag;
            end
        end
        if ~mod(j, 100)
            disp("Timestep: "+j+", time: "+toc)
        end
    end
    disp("Maximum tracking performance test completed, time:"+toc)
    completeness1fold = 100*sum(sum(visibilityMat>0))/(kTrg*size(tspan, 2));
    completeness3fold = 100*sum(sum(visibilityMat>2))/(kTrg*size(tspan, 2));
    disp("Completeness, 1-fold: "+completeness1fold+"%")
    disp("Completeness, 3-fold: "+completeness3fold+"%")
    
    %% 6) Number of observations
    observedMat = (visibilityMat > 0);
    observedMatHalfOrbitalPeriod = (visibilityMat > 0);
    for kTrg = 1:nTrg
        for j = 2:size(tspan, 2)
            if observedMat(kTrg, j-1) > 0
                observedMat(kTrg, j) = 1;
                observedMatHalfOrbitalPeriod(kTrg, j) = 1;
            end
            pastObservations = find(visibilityMat(kTrg, 1:(j-1))>0);
            if ~isempty(pastObservations)
                lastObserved = pastObservations(end);
                if visibilityMat(kTrg, j) == 0 && (j - lastObserved) > timestepsBeforeTrackLoss
                    observedMatHalfOrbitalPeriod(kTrg, j) = 0;
                end
            end
        end
    end
    nObservedTargets(kStrategy, :) = sum(observedMat);
    nObservedTargetsHalfOrbitalPeriod(kStrategy, :) = sum(observedMatHalfOrbitalPeriod);
    disp("Observation data processed, time: "+toc)
    
    %% 7) Track quality
    trackLengthVec = [];
    timeBetweenTracksVec = [];
    trackedTargetData = [];
    minObservations = size(tspan, 2);
    worstTarget = 0;
    for kTrg = 1:nTrg
        beginObservation = 0;
        endObservation = 0;
        nObservationsTemp = sum(visibilityMat(kTrg, :) > 0);
        if nObservationsTemp < minObservations
            worstTarget = kTrg;
            minObservations = nObservationsTemp;
        end
        for j = 1:size(tspan, 2)
            if j > 1 && visibilityMat(kTrg, j-1) > 0 && visibilityMat(kTrg, j) == 0
                endObservation = j;
                trackLength = (endObservation - beginObservation);
                if trackLength > minTrackLength
                    trackLengthVec = [trackLengthVec; trackLength * deltaT];
                    trackedTargetData = [trackedTargetData; [Target(kTrg).kep(j, :) Target(kTrg).size]];
                end
                beginObservation = 0;
            end
            if j == size(tspan, 2)
                endObservation = j;
                trackLength = (endObservation - beginObservation);
                if trackLength > minTrackLength
                    if beginObservation > 0
                        trackLengthVec = [trackLengthVec; trackLength * deltaT];
                    end
                    trackedTargetData = [trackedTargetData; [Target(kTrg).kep(j, :) Target(kTrg).size]];
                end
                beginObservation = 0;
            end
            if visibilityMat(kTrg, j) > 0 && beginObservation == 0
                beginObservation = j;
                if endObservation > 0
                    timeBetweenTracks = (beginObservation - endObservation);
                    timeBetweenTracksVec = [timeBetweenTracksVec; timeBetweenTracks*deltaT];
                end
            end
        end
    end
    averageTrackLength = mean(trackLengthVec/60);
    medianTrackLength = median(trackLengthVec/60);
    averageTimeBetweenTracks = mean(timeBetweenTracksVec/60);
    medianTimeBetweenTracks = median(timeBetweenTracksVec/60);
    percentageUsefulTrackPairs  = 100 *...
                                sum((timeBetweenTracksVec/deltaT < timestepsBeforeTrackLoss) .* ...
                                timeBetweenTracksVec/deltaT > minTrackSeparation) / size(timeBetweenTracksVec, 1);
    numUsefulTrackPairs         = sum((timeBetweenTracksVec/deltaT < timestepsBeforeTrackLoss) .* ...
                                timeBetweenTracksVec/deltaT > minTrackSeparation);
    numberOfTracks = size(trackLengthVec, 1);
    disp("Track data processed, time: "+toc)
    disp("Average track length: "+averageTrackLength+" minutes")
    disp("Average time between tracks: "+averageTimeBetweenTracks+" minutes")
    disp("Median track length: "+medianTrackLength+" minutes")
    disp("Median time between tracks: "+medianTimeBetweenTracks+" minutes")
    disp("Percentage of useful consecutive track pairs: "+percentageUsefulTrackPairs+"%")
    disp("Number of useful track pairs: "+numUsefulTrackPairs)
    disp("")
    
    %% Plots - single strategy
    
    % figure
    % subplot(2,1,1)
    % hold on
    % for k = 1:size(Constellation, 2)
    %     plot(tspan/ConstellationParameters.T, rad2deg(Constellation(k).orientation(:, 1)))
    % end
    % xlabel("time [T]")
    % ylabel("elevation [deg]")
    % subplot(2,1,2)
    % hold on
    % for k = 1:size(Constellation, 2)
    %     plot(tspan/ConstellationParameters.T, rad2deg(Constellation(k).orientation(:, 2)))
    % end
    % xlabel("time [T]")
    % ylabel("azimuth [deg]")

    % The following should be weighted against the original distributions!!!
    % figure()
    % subplot(2,2,1)
    % histogram(trackedTargetData(:, 1) - Constants.R_E)
    % xlabel("track sma height [km]")
    % subplot(2,2,2)
    % histogram(rad2deg(trackedTargetData(:, 3)))
    % xlabel("track inclination [deg]")
    % subplot(2,2,3)
    % histogram(rad2deg(trackedTargetData(:, 4)))
    % xlabel("track raan [deg]")
    % subplot(2,2,4)
    % histogram(rad2deg(log10(trackedTargetData(:, 7))))
    % xlabel("track size [m]")
    
    % figure()
    % subplot(1,2,1)
    % histogram(trackLengthVec/60, nBins)
    % xlabel("track length [min]")
    % subplot(1,2,2)
    % histogram(log10(trackLengthVec/60))
    % xlabel("log track length [lmin]")
    % 
    % figure()
    % subplot(1,2,1)
    % histogram(timeBetweenTracksVec/60, nBins)
    % xlabel("time between tracks [min]")
    % subplot(1,2,2)
    % histogram(log10(timeBetweenTracksVec/60))
    % xlabel("log time between tracks [lmin]")
    % 
    % figure()
    % imagesc(visibilityMat)
    % xlabel("timestep")
    % ylabel("target")
    % colorbar
    
    timestep = randi(size(tspan, 2));
    % figure()
    % axis equal
    % hold on
    % for k = 1:nObs
    %     scatter3(Constellation(k).x(timestep, 1),Constellation(k).x(timestep, 2),Constellation(k).x(timestep, 3), "^", "m", "filled")
    %     plotConicalFOV(VisibilityParams.halfFOV, 1e4, 50, Constellation(k).orientation(timestep,:)', Constellation(k).x(timestep,1:3)', Constellation(k).x(timestep,4:6)')
    % end
    % xlabel("X [km]")
    % ylabel("Y [km]")
    % zlabel("Z [km]")
    % title("Start")
end

%% Plots - strategy comparison
figure()
hold on
for k = 1:size(trackingStrategy,2)
    stairs(tspan/ConstellationParameters.T, 100*nObservedTargetsHalfOrbitalPeriod(k, :)/nTrg, "LineWidth", 2)
end
xlabel("time [T]", "FontSize", 12)
ylabel("percentage of targets obs. within T/2 [%]", "FontSize", 12)
l = legend(trackingStrategy);
l.Location = "southeast";

% 
% figure()
% hold on
% for kStrategy = 1:size(trackingStrategy, 2)
%     stairs(tspan/ConstellationParameters.T, 100*nObservedTargets(kStrategy, :)/nTrg)
% end
% xlabel("time [T]")
% ylabel("percentage of obs. targets [%]")
% l = legend(trackingStrategy);
% l.Location = "southeast";
% 
% figure()
% hold on
% for k = 1:size(trackingStrategy, 2)
%     stairs(tspan / ConstellationParameters.T, numberOfPotentialMeasurements(:, k))
% end
% xlabel("time [T]")
% ylabel("number of measurements")
% l = legend(trackingStrategy);
% l.Location = "southeast";

% figure()
% testTarget = Target(randi(nTrg));
% plot3(testTarget.x(:, 1),testTarget.x(:, 2),testTarget.x(:, 3))
% axis equal

if ~loadTargets
    figure()
    subplot(2,2,1)
    histogram(smaVec)
    xlabel("semi-major axis [km]")
    subplot(2,2,2)
    histogram(inclinationVec)
    xlabel("inclination [deg]")
    subplot(2,2,3)
    histogram(eccentricityVec)
    xlabel("eccentricity [-]")
    subplot(2,2,4)
    histogram(log10(eccentricityVec))
    xlabel("log eccentricity [-]")
end

% figure()
% scatter3(cart0Trg(:, 1), cart0Trg(:, 2), cart0Trg(:, 3), "filled")
% axis equal

disp(" ")
disp("Plot completed, time: "+toc)

%% FOV functions
function orientation    = getOrientationGreedy(prevOrientation, kT, Target, rObs, ...
                        DCM_ECI2LVLH, dirSun, R, sizeTrg, VisibilityParams, ...
                        maxSlewRate, nRuns)
    maxNumberOfMeasurements = 0;
    optimalOrientation = prevOrientation + (2*rand(2,1)-1)*maxSlewRate;
    for k = 1:nRuns
        numberOfMeasurements = 0;
        if k == 1
            currentOrientation = prevOrientation;
        else
            currentOrientation = prevOrientation + (2*rand(2,1)-1)*maxSlewRate;
        end
        for kTrg = 1:size(Target, 2)
            rTrg = Target(kTrg).x(kT, 1:3)';
            flag    = isTargetVisibleToObserver(rTrg, rObs, ...
                    getDirOrientation(currentOrientation), DCM_ECI2LVLH, dirSun, R, ...
                    sizeTrg, VisibilityParams);
            numberOfMeasurements = numberOfMeasurements + flag;
        end
        if numberOfMeasurements > maxNumberOfMeasurements
            maxNumberOfMeasurements = numberOfMeasurements;
            optimalOrientation = currentOrientation;
        end
    end
    orientation = optimalOrientation;
end

function orientation = getOrientation(strategy, prevOrientation, kT, Target, rObs, ...
                        DCM_ECI2LVLH, dirSun, R, sizeTrg, VisibilityParams, ...
                        maxSlewRate, nRuns)
    if nargin < 2 && strcmpi(strategy, "Greedy")
        strategy = "Along-track";
    end
    if nargin == 1 && strcmpi(strategy, "Greedy")
        strategy = "Along-track";
    end
    if strcmpi(strategy, "Radial")
        orientation = deg2rad([90; 0]);
    elseif strcmpi(strategy, "Along-track")
        orientation = deg2rad([90; 90]);
    elseif strcmpi(strategy, "Cross-track")
        orientation = deg2rad([0; 0]);
    elseif strcmpi(strategy, "Greedy")
        orientation     = getOrientationGreedy(prevOrientation, kT, Target, rObs, ...
                        DCM_ECI2LVLH, dirSun, R, sizeTrg, VisibilityParams, ...
                        maxSlewRate, nRuns);
    else
        orientation = [0; 0];
    end
end

function flag   = isWithinConicalFOV(observerPos, targetPos, dirOrientation, ...
                DCM_ECI2LVLH, halfFOV)
    dr_ECI = targetPos - observerPos;
    dr_LVLH = DCM_ECI2LVLH*dr_ECI;
    flag = ( acos(dot(dr_LVLH,dirOrientation)/norm(dr_LVLH)) <= halfFOV );
end

function dirOrientation = getDirOrientation(orientation)
    dirOrientation = nan(3,1);
    elv = orientation(1);
    az = orientation(2);
    dirOrientation(1,1) = sin(elv)*cos(az);
    dirOrientation(2,1) = sin(elv)*sin(az);
    dirOrientation(3,1) = cos(elv);
end

function DCM = get_DCM_ECI2LVLH(y_frame)
    % requires column vector
    r_vect = y_frame(1:3,1);
    v_vect = y_frame(4:6,1);
    h_vect = cross(r_vect,v_vect);
    r = norm(r_vect);
    h = norm(h_vect);
    e_x_lvlh = r_vect/r;
    e_z_lvlh = h_vect/h;
    e_y_lvlh = cross(e_z_lvlh, e_x_lvlh);
    DCM = [e_x_lvlh e_y_lvlh e_z_lvlh]';
end

function DCM = get_DCM_LVLH2LOCAL(elv, az)
    DCM = zeros(3,3);
    DCM(3,1) = sin(elv)*cos(az);
    DCM(3,2) = sin(elv)*sin(az);
    DCM(3,3) = cos(elv);
    DCM(1,1) = cos(elv)*cos(az);
    DCM(1,2) = cos(elv)*sin(az);
    DCM(1,3) = -sin(elv);
    DCM(2,1) = -sin(az);
    DCM(2,2) = cos(az);
end

function [] = plotConicalFOV(halfFOV, maxDist, N, orientation, observerPos, observerVel)
    % Make cone
    halfBase = tan(halfFOV);
    [X,Y,Z] = cylinder([0 halfBase], N);
    X = maxDist*X;
    Y = maxDist*Y;
    Z = maxDist*Z;
    
    % Get rotation matrix from local orientation to ECI
    DCM_ECI2LVLH = get_DCM_ECI2LVLH([observerPos; observerVel]);
    DCM_LVLH2LOCAL = get_DCM_LVLH2LOCAL(orientation(1),orientation(2));
    DCM_LOCAL2ECI = (DCM_LVLH2LOCAL*DCM_ECI2LVLH)';
    R = DCM_LOCAL2ECI;

    % Apply rotation + translation to observer position
    for k = 1:(N+1)
        X(1,k) = observerPos(1);
        Y(1,k) = observerPos(2);
        Z(1,k) = observerPos(3);
        temp = observerPos + R*[X(2,k);Y(2,k);Z(2,k)];
        X(2,k) = temp(1);
        Y(2,k) = temp(2);
        Z(2,k) = temp(3);
    end

    % Plot
    h=surf(X,Y,Z,'LineStyle','none','FaceAlpha',0.4);
end

%% Base tracking functions

function Target = propagateJ2(deltaT, Target)
    kepNext = Target.kep(end, :);
    e = kepNext(2);
    kepNext(4) = kepNext(4) + Target.raanDot * deltaT;
    kepNext(5) = kepNext(5) + Target.aopDot * deltaT;
    M = convert_anomaly_v2M(kepNext(6), e);
    M = M + Target.n * deltaT;
    kepNext(6) = convert_anomaly_M2v(M, e);
    Target.kep(end+1, :) = kepNext;
end

% Define visibility function
function flag = isTargetVisibleToObserver(rTrg, rObs, dirOrientation, DCM_ECI2LVLH, dirSun, R, trgSize, VisibilityParams)
    flag    = isAboveTheHorizon(rObs, rTrg, R) ...
            && ~isEclipsed(rTrg, dirSun, R) ...
            && ~isScatteringWeakMag(rObs, rTrg, dirSun, trgSize, VisibilityParams) ...
            && isWithinConicalFOV(rObs, rTrg, dirOrientation, DCM_ECI2LVLH, VisibilityParams.halfFOV);
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
