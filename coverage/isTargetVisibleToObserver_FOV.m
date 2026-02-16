function bool = isTargetVisibleToObserver_FOV(rTrg, rObs, dirSun, dirPointing, D_t, halfFOV, V_lim, R_H)
    if nargin < 8
        R_H = 6471;
    end
    if nargin < 7
        V_lim = 16;
    end
    if nargin < 6
        halfFOV = 0.06894; % FOV is 7.9 deg
    end
    bool = 0;
    if ~isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, V_lim, R_H)
        return
    end
    if ~isWithinFOV(rTrg, rObs, dirPointing, halfFOV)
        return
    end
    bool = 1;
end