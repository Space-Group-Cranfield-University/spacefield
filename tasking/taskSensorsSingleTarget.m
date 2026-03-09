function OBS = taskSensorsSingleTarget...
            (timestep, OBS, TRG)
    if size(TRG, 2) > 1
        error("Single target tasking strategy requires a single target. Select target or change tasking options")
    end
    rTrg = TRG.xMat(timestep, 1:3)';
    for k = 1:size(OBS, 2)
        rObs = OBS(k).xMat(timestep, 1:3)';
        dr = rTrg - rObs;
        prePoint = OBS(k).pointingMat(timestep - 1, :)';
        newPoint_ECI = dr / norm(dr);
        DCM_LOCAL2ECI = getRotationMatrixLocalToECI(OBS(k).xMat(timestep, :)');
        newPoint = DCM_LOCAL2ECI' * newPoint_ECI;
        requiredSlewAngle = acos(prePoint' * newPoint);
        if requiredSlewAngle > OBS(k).SensorParameters.maxSlewAngle
            rotationAxis = cross(prePoint, newPoint);
            theta = OBS(k).SensorParameters.maxSlewAngle;
            newPoint = applyRodriguezFormula(prePoint, rotationAxis, theta);
        end
        raDec = convertPositionToRaDec(newPoint);
        %if raDec(end) > OBS(k).SensorParameters.maxEl
        %    newPoint = prePoint;
        %end
        OBS(k).pointingMat(timestep, :) = newPoint';
    end
end