function OPTIONS = getStandardSimulationOptions()
    OPTIONS.verbose = 1;
    OPTIONS.taskingStrategy = "static";
    OPTIONS.filterFlag = 0;
    OPTIONS.postProcessing = 1;
    OPTIONS.getCoverage = 0;
    OPTIONS.maxOrbitalFraction = 0.5;
    OPTIONS.minTrackLength = 120;
end