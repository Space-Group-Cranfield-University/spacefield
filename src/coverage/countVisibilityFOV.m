function n = countVisibilityFOV(rTrg, rObsMat, dirPointingMat, dirSun, D_t, SensorParameters, R_H)
    if nargin < 6
        R_H = 6471;
    end
    if nargin < 5
        SensorParameters = getReducedSensorParameters;
    end
    n = 0;
    for k = 1:size(rObsMat, 1)
        rObs = rObsMat(k, :)';
        dirPointing = dirPointingMat(k, :)';
        if isTargetVisibleToObserver_FOV(rTrg, rObs, dirPointing, dirSun, D_t, SensorParameters, R_H)
            nVisible = nVisible + 1;
        end
    end
end