function OPTIONS = getStandardRadarSearchOptions()
    OPTIONS.nRuns = 100;
    OPTIONS.RadarParameters = getStandardRadar;
    OPTIONS.dirSun = [-1 0 0]';
    OPTIONS.nFold = 1;
    OPTIONS.nTrg = 1e4;
end