function AzEl = getOpticalMeasurement(xTrg, xObs)
    rTrg = xTrg(1:3);
    rObs = xObs(1:3);
    deltaR = rTrg - rObs;
    DCM = getRotationMatrixLVLHToECI(xObs);
    deltaR = DCM' * deltaR;
    Az = atan2(deltaR(2), deltaR(1));
    El = asin(deltaR(3) / norm(deltaR));
    AzEl = [Az; El];
end