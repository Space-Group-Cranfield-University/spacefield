function OPTIONS = getStandardConstellationSearchOptions()
    OPTIONS.dN = 0;
    OPTIONS.nRuns = 100;
    OPTIONS.nTests = 1;
    OPTIONS.nFold = 1;
    OPTIONS.dirSun0 = [-1 0 0]';
    OPTIONS.n_t_season = 4;
    OPTIONS.n_t_day = 3;
    OPTIONS.isSeasonal = 0;
    OPTIONS.isDaily = 0;
    OPTIONS.nTrg = 1e4;
    OPTIONS.SensorParameters = getReducedSensorParameters;
end