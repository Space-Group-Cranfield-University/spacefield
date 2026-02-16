function bool = isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, V_lim, R_H)
    if nargin < 6
        R_H = 6471;
    end
    if nargin < 5
        V_lim = 16;
    end
    bool = 0;
    if isEclipsed(rTrg, dirSun, R_H)
        return
    end
    if ~isAboveTheHorizon(rObs, rTrg, R_H)
        return
    end
    if isScatteringWeak(rObs, rTrg, dirSun, D_t, V_lim)
        return
    end
    bool = 1;
end