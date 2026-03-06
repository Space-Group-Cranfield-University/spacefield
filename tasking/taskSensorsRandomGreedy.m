function OBS = taskSensorsRandomGreedy(timestep, OBS, TRG, prevObsTrgCrossVisibilityMat, dirSun, OPTIONS)
    % Requires pointing directions to be initialized within the bounds
    % given by OPTIONS.
    if nargin < 5 || isempty(OPTIONS)
        OPTIONS = getStandardGreedyOptions;
    end
    
    if OPTIONS.optimizeForCompleteNetwork
        OBS = taskSensorsRandomGreedyNetwork...
        (timestep, OBS, TRG, prevObsTrgCrossVisibilityMat, dirSun, OPTIONS);
    else
        OBS = taskSensorsRandomGreedySingleSensor...
        (timestep, OBS, TRG, prevObsTrgCrossVisibilityMat, dirSun, OPTIONS);
    end
end