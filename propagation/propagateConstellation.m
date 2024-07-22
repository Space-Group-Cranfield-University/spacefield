% propagateConstellation propagates a given Walker constellation according
% to Keplerian dynamics.
%
%   propagateConstellation(timeVec, Constellation, Constants)
%
%   Inputs:
%   - timeVec       : row vector containing the propagation timesteps [s].
%   - Constellation : constellation structure (see InitialiseConstellation.m)
%   - Constants     : astronomical constants (see InitialiseAstronomicalConstants.m)
%
%   Outputs:
%   - Constellation : updated constellation structure.

function Constellation = propagateConstellation(timeVec, Constellation, Constants)
    T = computeOrbitalPeriod(Constellation(1).a, Constants.MU_E);
    thetaVec = 2*pi * timeVec / T;
    for iTheta = 1:size(thetaVec, 2)
        for iSat = 1:size(Constellation, 2)
            currentKep = [Constellation(iSat).kep(1:5, 1); Constellation(iSat).kep(6, 1) + thetaVec(iTheta)];
            Constellation(iSat).x(iTheta, :) = convertKepToCart(currentKep, Constants.MU_E)';
        end
    end
end