% convertRadToArcsec Converts radians to arcseconds.
%
%   angleArcsec = convertRadToArcsec(angleRad) converts an angle angleRad
%   in radians to an angle angleArcsec in arcseconds.

function angleArcsec = convertRadToArcsec(angleRad)
    angleArcsec = angleRad * 206265;
end