function bool = isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, SensorParameters, R_H)
    if nargin < 6
        R_H = 6471;
    end
    if nargin < 5
        SensorParameters = getReducedSensorParameters();
    end
    bool = 0;
    if isEclipsed(rTrg, dirSun, R_H)
        return
    end
    if ~isAboveTheHorizon(rObs, rTrg, SensorParameters.alpha_e, R_H)
        return
    end
    if ~isSunExcluded(rObs, rTrg, dirSun, SensorParameters.alpha_s)
        return
    end
    if isScatteringWeak(rObs, rTrg, dirSun, D_t, SensorParameters.V_lim)
        return
    end
    bool = 1;
end