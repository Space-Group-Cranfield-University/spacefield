% isObstructed checks wheter the line of sight connecting an observer to a
% target object is obstructed by a central body.
%
%   flag = isObstructed(rObs, rTrg, R)
%
%   Inputs:
%   - rObs      : displacement vector of the observer.
%   - rTrg      : displacement vector of the target object.
%   - R         : radius of central body plus atmospheric height.
%
%   Outputs:
%   - flag      : true if the line of sight is obstructed, false otherwise.


function flag = isAboveTheHorizon(rObs, rTrg, R)

    dr = rObs - rTrg;
    dr_norm = norm(dr);
    r_obs_surf = sqrt(rObs' * rObs - R^2);
    relative_visibility = dr_norm * r_obs_surf - dr' * rObs;
    flag = relative_visibility > 0;

end