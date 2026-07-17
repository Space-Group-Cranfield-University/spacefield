function OBS = taskSensors...
    (timestep, OBS, taskingStrategy, TRG, prevObsTrgCrossVisibilityMat, ...
    dirSun, OPTIONS)
    if nargin < 7
        OPTIONS = [];
    end
    if nargin < 6
        dirSun = [];
    end
    if nargin < 5
        prevObsTrgCrossVisibilityMat = [];
    end
    if nargin < 4
        TRG = [];
    end
    if nargin < 3
        taskingStrategy = "static";
    end

    if strcmp(taskingStrategy, "static")
        OBS = taskSensorsStatic(timestep, OBS);
    elseif strcmp(taskingStrategy, "random-greedy")
        OBS = taskSensorsRandomGreedy...
            (timestep, OBS, TRG, prevObsTrgCrossVisibilityMat, dirSun, OPTIONS);
    elseif strcmp(taskingStrategy, "single-target")
        OBS = taskSensorsSingleTarget...
            (timestep, OBS, TRG);
    end
end