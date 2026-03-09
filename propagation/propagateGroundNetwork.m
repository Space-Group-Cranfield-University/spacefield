function GND = propagateGroundNetwork(timeVec, GND)
    OM = initializeAstronomicalConstants().OM_E;
    omVec = [0; 0; OM];
    for k = 1:size(GND, 2)
        for j = 1:size(timeVec, 2)
            DCM = getRotationMatrixECEFToECI(timeVec(j) - GND(k).t_0);
            r_ECI = DCM * GND(k).r0_ECEF;
            v_ECI = cross(omVec, r_ECI);
            GND(k).xMat(j, :) = [r_ECI', v_ECI'];
        end
    end
end