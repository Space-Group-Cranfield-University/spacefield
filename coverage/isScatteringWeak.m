% isBlinded checks whether an observer is blinded by the Sun, i.e. whether
% the target is not reflecting sunlight back to the observer due to their
% relative geometry. If the observer is behind the target with respect to
% the Sun, it cannot receive the light that is reflected from the target,
% and therefore it cannot detect the target.
%
%   flag = isBlinded(rObs, rTrg, dirSun)
%
%   Inputs:
%   - rObs      : displacement vector of the observer.
%   - rTrg      : displacement vector of the target object.
%   - dirSun    : unit vector pointing towards the Sun.
%
%   Outputs:
%   - flag      : true if the observer is blinded, false otherwise.

function flag = isScatteringWeak(rObs, rTrg, dirSun)

    dr = rObs - rTrg;
    sun_reflection = dr' * dirSun;
    flag = sun_reflection < 0;

end