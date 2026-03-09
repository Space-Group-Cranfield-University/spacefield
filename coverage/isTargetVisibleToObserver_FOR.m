function bool = isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, SensorParameters)
    if nargin < 5
        SensorParameters = getReducedSensorParameters();
    end
    bool = 0;
    if SensorParameters.sensorType == "optical" && isEclipsed(rTrg, dirSun, SensorParameters.R_H)
        return
    end
    if ~isAboveTheHorizon(rObs, rTrg, SensorParameters.alpha_e, SensorParameters.R_H)
        return
    end
    if SensorParameters.sensorType == "optical" && ~isSunExcluded(rObs, rTrg, dirSun, SensorParameters.alpha_s)
        return
    end
    if isScatteringWeak(rObs, rTrg, dirSun, D_t, SensorParameters)
        return
    end
    bool = 1;
end