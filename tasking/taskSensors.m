function OBS = taskSensors(timestep, OBS, taskingStrategy)
    if nargin < 3
        taskingStrategy = "static";
    end

    if strcmp(taskingStrategy, "static")
        for k = 1:size(OBS, 2)
            OBS(k).pointingMat(timestep, :) = OBS(k).pointingMat(timestep - 1, :);
        end
    end
end