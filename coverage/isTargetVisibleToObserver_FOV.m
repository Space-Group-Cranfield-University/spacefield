function bool = isTargetVisibleToObserver_FOV(rTrg, rObs, dirPointing, dirSun, D_t, SensorParameters)
    if nargin < 6
        SensorParameters = getReducedSensorParameters;
    end
    bool = 0;
    if ~isfield(SensorParameters, "pointingDirBody")
        for k = 1:size(SensorParameters, 2)
            SensorParameters(k).pointingDirBody = [0 0]';
        end
    end
    if ~isWithinFOV...
            (rTrg, rObs, ...
            getCompletePointingDirection(dirPointing, SensorParameters.pointingDirBody), ...
            SensorParameters.halfFOV)
        return
    end
    if ~isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, SensorParameters)
        return
    end
    bool = 1;
end