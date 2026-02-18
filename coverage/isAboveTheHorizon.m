% isObstructed checks wheter the line of sight connecting an observer to a
% target object is obstructed by a central body.
%
%   flag = isObstructed(rObs, rTrg, R)
%
%   Inputs:
%   - rObs      : displacement vector of the observer.
%   - rTrg      : displacement vector of the target object.
%   - alpha_e   : Earth's exclusion angle
%   - R         : radius of central body plus atmospheric height.
%
%   Outputs:
%   - flag      : true if the line of sight is obstructed, false otherwise.


function flag = isAboveTheHorizon(rObs, rTrg, alpha_e, R)
    if nargin < 4
        R = 6371;
    end
    if nargin < 3
        alpha_e = deg2rad(25);
    end
    rObsNorm = norm(rObs);
    theta_e = arcsin( R / rObsNorm );
    dr = rObs - rTrg;
    dr_norm = norm(dr);
    target_angle = arccos( - rObs' * dr / (rObsNorm * dr_norm) );
    relative_visibility = target_angle - theta_e - alpha_e;
    flag = (relative_visibility > 0);

end