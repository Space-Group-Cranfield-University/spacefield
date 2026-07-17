function DCM = getRotationMatrixLocalToECI(x)
    DCM_LocalToLVLH = getRotationMatrixLocalToLVLH;
    DCM_LVLHToECI = getRotationMatrixLVLHToECI(x);
    DCM = DCM_LVLHToECI * DCM_LocalToLVLH;
end