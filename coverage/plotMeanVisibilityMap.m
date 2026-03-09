function f = plotMeanVisibilityMap...
    (timeVec, OBS, dirSun, h_trg, D_trg, SensorParameters, ...
    plotGroundTracksFlag, N_radec, n_t_ground_track)
    if nargin < 9
        n_t_ground_track = 50;
    end
    if nargin < 8
        N_radec = 2*1e2;
    end
    if nargin < 7
        plotGroundTracksFlag = 0;
    end
    if nargin < 2
        OBS = initializeWalkerConstellation;
    end
    if nargin < 6
        SensorParameters = getReducedSensorParameters;
    end
    if nargin < 5
        D_trg = 1;
    end
    if nargin < 4
        h_trg = 800;
    end
    if nargin < 3
        dirSun = [-1 0 0]';
    end
    meanVisibilityMat = 0;
    if ~isfield(OBS, "GND")
        if ~isfield(timeVec, "Day")
            for j = 1:size(timeVec, 2)
                rObsMat = getConstellationPositionMatrix(OBS, j);
                dirSunTemp = dirSun;
                [raMat, decMat, visibilityCountMat] = getVisibilityGrid(rObsMat, dirSunTemp, SensorParameters, D_trg, h_trg, N_radec);
                meanVisibilityMat = meanVisibilityMat + visibilityCountMat;
            end
            N_t = size(timeVec, 2);
        else
            dirSunMat = propagateDirSun(timeVec.Season, dirSun);
            for k = 1:size(timeVec.Season, 2)
                timeVecTemp = timeVec.Season(k) + timeVec.Day;
                OBS = propagateConstellation(timeVecTemp, OBS);
                dirSunTemp = dirSunMat(k, :)';
                for j = 1:size(timeVec.Day, 2)
                    rObsMat = getConstellationPositionMatrix(OBS, j);
                    [raMat, decMat, visibilityCountMat] = ...
                        getVisibilityGrid...
                        (rObsMat, dirSunTemp, SensorParameters, D_trg, h_trg, N_radec);
                    meanVisibilityMat = meanVisibilityMat + visibilityCountMat;
                end
                disp("Day finished: "+string(k)+" / "+string(size(timeVec.Season, 2)))
            end
            N_t = size(timeVec.Day, 2) * size(timeVec.Season, 2);
        end
    else
        for k = 1:size(OBS, 2)
            currentGND = OBS(k).GND;
            rObsMat = getConstellationPositionMatrix(currentGND);
            dirSunTemp = dirSun;
            [raMat, decMat, visibilityCountMat] = ...
                getVisibilityGrid...
                (rObsMat, dirSunTemp, SensorParameters, D_trg, h_trg, N_radec);
            meanVisibilityMat = meanVisibilityMat + visibilityCountMat;
        end
        N_t = size(OBS, 2);
    end
    meanVisibilityMat = meanVisibilityMat / N_t;
    if SensorParameters.sensorType == "optical"
        timeVec_gt = linspace(0, OBS(1).T, n_t_ground_track);
        if ~isfield(OBS, "raDecMat")
            OBS = propagateConstellation(timeVec_gt, OBS);
            OBS = getConstellationRaDec(OBS);
        end
        f = figure();
        m = mesh(rad2deg(raMat), rad2deg(decMat), meanVisibilityMat);
        hold on
        if plotGroundTracksFlag
            for k = 1:size(OBS, 2)
                plotGroundTracks(OBS(k).raDecMat)
            end
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
    else
        f = figure();
        m = mesh(rad2deg(raMat), rad2deg(decMat), meanVisibilityMat);
        hold on
        if plotGroundTracksFlag
            OBS = getConstellationRaDec(OBS);
            for k = 1:size(OBS, 2)
                plotGroundTracks(OBS(k).raDecMat)
            end
        end
        m.FaceColor = 'flat';
        view(2)
        axis equal
        xlabel("Right Ascension [deg]", "FontSize", 14);
        ylabel("Declination [deg]", "FontSize", 14);
        title("Visibility Map, h_t = "+string(h_trg)+" km, D_t = "+string(D_trg)+" m", "FontSize", 12);
        colorbar;
        set(gca, "FontSize", 14)
    end
end