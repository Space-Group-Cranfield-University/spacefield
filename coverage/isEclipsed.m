% isEclipsed checks whether a given target object is shadowed by a
% central body (for instance, the Earth), i.e. is in eclipse condition.
%
%   flag = isEclipsed(rTrg, dirSun, R)
%
%   Inputs:
%   - rTrg      : displacement vector of the target object.
%   - dirSun    : unit vector pointing towards the Sun.
%   - R         : radius of central body plus atmospheric height.
%
%   Outputs:
%   - flag      : true if the target is eclipsed, false otherwise.


function flag = isEclipsed(rTrg, dirSun, R)
    if nargin < 3
        R = 6471;
    end

    % positive if the target is behind the plane of the terminator
    eclipse_behind_earth = - rTrg' * dirSun;

    % positive if the target is within the cylinder tangent to the Earth 
    % at the terminator
    eclipse_within_earth = (rTrg' * dirSun)^2 - norm(rTrg)^2 + R^2;

    flag = (eclipse_behind_earth > 0) && (eclipse_within_earth > 0);

end