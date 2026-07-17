function bool = isTargetVisibleToConstellation_FOV(rTrg, rObsMat, dirPointingMat, dirSun, D_t, nFold, SensorParameters, R_H)
    if nargin < 8
        R_H = 6471;
    end
    if nargin < 7
        SensorParameters = getReducedSensorParameters;
    end
    if nargin < 6
        nFold = 1;
    end
    bool = 0;
    nVisible = 0;
    for k = 1:size(rObsMat, 1)
        rObs = rObsMat(k, :)';
        dirPointing = dirPointingMat(k, :)';
        if isTargetVisibleToObserver_FOV(rTrg, rObs, dirPointing, dirSun, D_t, SensorParameters, R_H)
            nVisible = nVisible + 1;
        end
        if nVisible == nFold
            bool = 1;
            break
        end
    end
end