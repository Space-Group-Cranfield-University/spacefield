function f = plotVisibilityMap(OBS, dirSun, h_trg, D_trg, SensorParameters, rObsMat, N_radec, n_t)
    if nargin < 8
        n_t = 50;
    end
    if nargin < 7
        N_radec = 2*1e2;
    end
    if nargin < 1
        OBS = initializeWalkerConstellation;
    end
    if nargin < 6
        rObsMat = getConstellationPositionMatrix(OBS);
    end
    if nargin < 5
        SensorParameters = getReducedSensorParameters;
    end
    if nargin < 4
        D_trg = 1;
    end
    if nargin < 3
        h_trg = 800;
    end
    if nargin < 2
        dirSun = [-1 0 0]';
    end
    timeVec = linspace(0, OBS(1).T, n_t);
    [raMat, decMat, visibilityCountMat] = getVisibilityGrid(rObsMat, dirSun, SensorParameters, D_trg, h_trg, N_radec);
    OBS = propagateConstellation(timeVec, OBS);
    OBS = getConstellationRaDec(OBS);
    f = figure();
    m = mesh(rad2deg(raMat), rad2deg(decMat), visibilityCountMat);
    hold on
    for k = 1:size(OBS, 2)
        plotGroundTracks(OBS(k).raDecMat)
    end
    m.FaceColor = 'flat';
    view(2)
    axis equal
    xlabel("Right Ascension [deg]", "FontSize", 14);
    ylabel("Declination [deg]", "FontSize", 14);
    %zlabel('Visibility Count');
    title("Visibility Map, h_t = "+string(h_trg)+" km, D_t = "+string(D_trg)+" m", "FontSize", 12);
    colorbar;
    set(gca, "FontSize", 14)
end