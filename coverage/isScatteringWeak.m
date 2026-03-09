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

function flag = isScatteringWeak(rObs, rTrg, dirSun, A_t, SensorParameters, isCrossSection)
    if nargin < 6
        isCrossSection = 0;
    end
    if ~isCrossSection
        D_t = A_t;
        A_t = A_t^2;
    else
        D_t = sqrt(A_t);
    end
    dr = rObs - rTrg;
    deltaR = norm( dr );
    alpha = acos( dr' * dirSun / deltaR );
    if SensorParameters.sensorType == "optical"
        V = computeApparentMagnitude( deltaR*1e3, A_t, alpha );
        flag = (V > SensorParameters.V_lim);
    elseif SensorParameters.sensorType == "radar"
        SNR = computeRadarSNR(deltaR * 1e3, D_t, SensorParameters);
        flag = (SNR < SensorParameters.minSNR);
    end
end