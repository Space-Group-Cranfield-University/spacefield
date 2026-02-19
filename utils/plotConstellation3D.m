function f = plotConstellation3D(OBS, rObsMat, orbitFlag, satFlag, az, el, n_t)
    f = figure();
    if nargin < 7
        n_t = 50;
    end
    if nargin < 6
        el = 30;
    end
    if nargin < 5
        az = 135;
    end
    if nargin < 4
        satFlag = 1;
    end
    if nargin < 3
        orbitFlag = 1;
    end
    if nargin < 1
        OBS = initializeWalkerConstellation;
    end
    if nargin < 2
        rObsMat = getConstellationPositionMatrix(OBS);
    end
    timeVec = linspace(0, OBS(1).T, n_t);
    OBS = propagateConstellation(timeVec, OBS);

    axis equal
    %set(gca,'Color','k')
    hold on
    xlabel("X, ECI [km]", "FontSize", 12)
    ylabel("Y, ECI [km]", "FontSize", 12)
    zlabel("Z, ECI [km]", "FontSize", 12)
    
    % Earth surface
    cdata = imread("2k_earth_daymap_clouds.jpg");
    [x,y,z] = sphere(100);
    %z = 0.001*z;
    Re = initializeAstronomicalConstants().R_E;
    earth = surf(Re*x, Re*y, -Re*z);
    set(earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none'); % "texturemap"
    %surf(X, Y, Z, "EdgeColor", "none", "FaceColor", "b")
    
    % Orbits
    if orbitFlag
        for iOrbits = 1:size(OBS, 2)
            plot3(OBS(iOrbits).xMat(:,1), OBS(iOrbits).xMat(:,2), OBS(iOrbits).xMat(:,3), "Color", "#EDB120", "LineWidth", 1)
        end
    end
    
    % Satellites
    if satFlag
        scatter3(rObsMat(:,1), rObsMat(:,2), rObsMat(:,3), "filled", "v", "MarkerFaceColor", "#D95319")
    end
    
    %lighting phong
    view(az, el)
    axis equal
    set(gca, "FontSize", 12)
    title("Constellation 3D plot", "FontSize", 12)
end