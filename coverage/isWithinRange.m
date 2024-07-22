% isWithinRange checks whether a target is within the visibility range of
% an observer.
%
%   flag = isWithinRange(rObs, rTrg, deltaR)
%
%   Inputs:
%   - rObs      : displacement vector of the observer.
%   - rTrg      : displacement vector of the target object.
%   - deltaR    : visibility range of sensor.
%
%   Outputs:
%   - flag      : true if the target is within range, false otherwise.


function flag = isWithinRange(rObs, rTrg, deltaR)

    flag = norm(rObs - rTrg) <= deltaR;

end