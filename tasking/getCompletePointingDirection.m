function dirPointingFull = getCompletePointingDirection(dirPointing, pointingDirBody)
    % Given a pointing direction of a sensor in the Body frame, and the
    % poitning direction of the spacecraft in the Local frame, outputs 
    % pointing direction of the sensor in Local frame.
    pointingVectorBody = convertRaDecToPosition([1, pointingDirBody']');
    DCM_BODY2LOCAL_Y = getRotationMatrix_Y(dirPointing(2))';
    DCM_BODY2LOCAL_Z = getRotationMatrix_Z(dirPointing(1));
    pointingVectorLocal = DCM_BODY2LOCAL_Z * DCM_BODY2LOCAL_Y * pointingVectorBody;
    dirPointingFull = convertPositionToRaDec(pointingVectorLocal);
    dirPointingFull = dirPointingFull(2:3);
end