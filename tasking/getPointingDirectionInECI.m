function dirPointing = getPointingDirectionInECI(currentObserver, timestep)
    if nargin < 2
        timestep = 1;
    end
    pointingDir = currentObserver.pointingMat(timestep, :)';
    if strcmp(currentObserver.pointingFrame, "Local")
        DCM_1 = getRotationMatrixLocalToLVLH;
        DCM_2 = getRotationMatrixLVLHToECI(currentObserver.xMat(timestep, :)');
        DCM = DCM_2 * DCM_1;
    elseif strcmp(currentObserver.pointingFrame, "LVLH")
        DCM = getRotationMatrixLVLHToECI(currentObserver.xMat(timestep, :)');
    else
        DCM = eye(3);
    end
    dirPointing = DCM * pointingDir;
end