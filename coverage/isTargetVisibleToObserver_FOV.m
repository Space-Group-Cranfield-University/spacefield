function bool = isTargetVisibleToObserver_FOV(rTrg, rObs, dirPointing, dirSun, D_t, SensorParameters)
    if nargin < 6
        SensorParameters = getReducedSensorParameters;
    end
    bool = 0;
    if ~isWithinFOV(rTrg, rObs, dirPointing, SensorParameters.halfFOV)
        return
    end
    if ~isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, SensorParameters)
        return
    end
    bool = 1;
end