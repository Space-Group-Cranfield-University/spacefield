function bool = isSunExcluded(rObs, rTrg, dirSun, alpha_s, R)
    if nargin < 5
        R = 6371;
    end
    if nargin < 4
        alpha_s = deg2rad(25);
    end
    bool = 1;
    if isEclipsed(rObs, dirSun, R)
        return
    end
    dr = rTrg - rObs;
    drNorm = norm(dr);
    targetAngle = arccos( dr'*dirSun / drNorm );
    sunExclusionDifference = targetAngle - alpha_s;
    bool = sunExclusionDifference > 0;
end