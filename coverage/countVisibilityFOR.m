function nVisible = countVisibilityFOR(rTrg, rObsMat, dirSun, D_t, SensorParameters)
    if nargin < 5
        SensorParameters = getReducedSensorParameters;
    end
    nVisible = 0;
    for k = 1:size(rObsMat, 1)
        rObs = rObsMat(k, :)';
        if isTargetVisibleToObserver_FOR(rTrg, rObs, dirSun, D_t, SensorParameters)
            nVisible = nVisible + 1;
        end
    end
end