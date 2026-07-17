function [GND, rObsMat] = initializeGroundNetwork(gndLocationMat, SensorParameters, t_0)
    if nargin < 3
        t_0 = 0;
    end
    R = initializeAstronomicalConstants().R_E;
    OM = initializeAstronomicalConstants().OM_E;
    om = [0, 0, OM]';
    for k = 1:size(gndLocationMat, 1)
        GND(k).t_0 = t_0;
        GND(k).lonLat = gndLocationMat(k, :);
        GND(k).SensorParameters = SensorParameters;
        GND(k).r0_sph = [R, GND(k).lonLat]';
        GND(k).r0_ECEF = convertRaDecToPosition(GND(k).r0_sph');
        DCM_ECEF2ECI = getRotationMatrixECEFToECI(t_0);
        GND(k).r0_ECI = DCM_ECEF2ECI * GND(k).r0_ECEF;
        GND(k).v0 = cross(om, GND(k).r0_ECI);
        GND(k).x0 = [GND(k).r0_ECI; GND(k).v0];
    end
    GND = setInitialPointingDirections(0, GND);
    rObsMat = getConstellationPositionMatrix(GND);
end