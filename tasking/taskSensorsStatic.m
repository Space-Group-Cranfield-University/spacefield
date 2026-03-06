function OBS = taskSensorsStatic(timestep, OBS)
    for k = 1:size(OBS, 2)
        OBS(k).pointingMat(timestep, :) = OBS(k).pointingMat(timestep - 1, :);
    end
end