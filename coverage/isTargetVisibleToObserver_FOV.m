function bool = isTargetVisibleToObserver_FOV(rTrg, rObs, dirPointing, dirSun, D_t, SensorParameters, R_H)
    if nargin < 7
        R_H = 6471;
    end
    if nargin < 6
        SensorParameters = getReducedSensorParameters;
    end
    bool = 0;
    if ~isWithinFOV(rTrg, rObs, dirPointing, SensorParameters.halfFOV)
        return
    end
    if ~isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, SensorParameters, R_H)
        return
    end
    bool = 1;
end