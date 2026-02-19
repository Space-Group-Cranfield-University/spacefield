function bool = isSunExcluded(rObs, rTrg, dirSun, alpha_s)
    if nargin < 4
        alpha_s = deg2rad(25);
    end
    bool = 1;
    if isEclipsed(rObs, dirSun)
        return
    end
    dr = rTrg - rObs;
    drNorm = norm(dr);
    targetAngle = acos( dr'*dirSun / drNorm );
    sunExclusionDifference = targetAngle - alpha_s;
    bool = (sunExclusionDifference > 0);
end