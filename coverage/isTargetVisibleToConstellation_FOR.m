function bool = isTargetVisibleToConstellation_FOR(rTrg, rObsMat, dirSun, D_t, nFold, V_lim, R_H)
    if nargin < 7
        R_H = 6471;
    end
    if nargin < 6
        V_lim = 16;
    end
    if nargin < 5
        nFold = 1;
    end
    bool = 0;
    nVisible = 0;
    for k = 1:size(rObsMat, 1)
        rObs = rObsMat(k, :)';
        if isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, V_lim, R_H)
            nVisible = nVisible + 1;
        end
        if nVisible == nFold
            bool = 1;
        end
    end
end