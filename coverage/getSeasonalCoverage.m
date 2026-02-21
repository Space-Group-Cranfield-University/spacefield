function [meanCoverage, varCoverage, seasonalCoverageVec] = getSeasonalCoverage(OBS, dirSun0, n_t_season, n_t_day_season, nFold)
    yearInSeconds = 365*86400;
    timeVecSeason = linspace(0, yearInSeconds, n_t_season + 1);
    timeVecSeason = timeVecSeason(1:(end-1));
    if nargin < 4
        n_t_day_season = 5;
    end
    timeVecDailyWithinSeason = linspace(0, OBS(1).T, n_t_day_season);
    timeVecDailyWithinSeason = timeVecDailyWithinSeason(1:(end-1));
    OBS_Season = propagateConstellation(timeVecSeason, OBS, 1);
    OBS_Day = OBS_Season;
    dirSunMat = propagateDirSun(timeVecSeason, dirSun0);
    seasonalCoverageVec = zeros(1, size(timeVecSeason, 2));
    for j = 1:size(timeVecSeason, 2)
        % Set initial day time to the j-th timestep of the season
        for k = 1:size(OBS_Day, 2)
            OBS_Day(k).x0 = OBS_Season(k).xMat(j, :)';
            OBS_Day(k).xMat = [];
        end
        OBS_Day = propagateConstellation(timeVecDailyWithinSeason, OBS_Day);
        dailyCoverageVec = getDailyCoverage(timeVecDailyWithinSeason, OBS_Day, dirSunMat(j, :)', nFold);
        seasonalCoverageVec(j) = mean(dailyCoverageVec);
    end
    meanCoverage = mean(seasonalCoverageVec);
    varCoverage = var(seasonalCoverageVec);
end