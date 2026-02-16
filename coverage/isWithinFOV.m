function bool = isWithinFOV(rTrg, rObs, dirPointing, halfFOV)
    dr = rTrg - rObs;
    theta = acos(dr' * dirPointing / norm(dr));
    bool = (theta < halfFOV)
end