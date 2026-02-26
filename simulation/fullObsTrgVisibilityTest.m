function ObsTrgCrossVisibilityMat = fullObsTrgVisibilityTest(timestep, OBS, TRG, dirSunMat, SensorParameters, R_H)
    if nargin < 6
        R_H = 100 + initializeAstronomicalConstants().R_E;
    end
    dirSun = dirSunMat(timestep, :)';
    ObsTrgCrossVisibilityMat = zeros(size(OBS, 2), size(TRG, 2));
    for j = 1:size(TRG, 2)
        rTrg = TRG(j).xMat(timestep, 1:3)';
        D_t = TRG(j).D_t;
        for k = 1:size(OBS, 2)
            rObs = OBS(k).xMat(timestep, 1:3)';
            dirPointing = getPointingDirectionInECI(OBS(k), timestep);
            ObsTrgCrossVisibilityMat(k, j) = isTargetVisibleToObserver_FOV(rTrg, rObs, dirPointing, dirSun, D_t, SensorParameters, R_H);
        end
    end
end