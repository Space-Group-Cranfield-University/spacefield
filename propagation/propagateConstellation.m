% propagateConstellation propagates a given Walker constellation according
% to Keplerian dynamics.
%
%   propagateConstellation(timeVec, Constellation, MU_E)
%
%   Inputs:
%   - timeVec       : row vector containing the propagation timesteps [s].
%   - Constellation : constellation structure (see InitialiseConstellation.m)
%   - MU_E          : gravitational parameter central body (see InitialiseAstronomicalConstants.m)
%
%   Outputs:
%   - Constellation : updated constellation structure.
%
%   NOTICE: ONLY works for circular constellation orbits

function Constellation = propagateConstellation(timeVec, Constellation, mu)
    if nargin < 3
        mu = initializeAstronomicalConstants().MU_E;
    end
    T = computeOrbitalPeriod(Constellation(1).a, mu);
    thetaVec = 2*pi * timeVec / T;
    for iSat = 1:size(Constellation, 2)
        Constellation(iSat).xMat = zeros(size(thetaVec, 2), 6);
        for iTheta = 1:size(thetaVec, 2) 
            currentKep = [Constellation(iSat).kep(1:5, 1); Constellation(iSat).kep(6, 1) + thetaVec(iTheta)];
            Constellation(iSat).xMat(iTheta, :) = convertKepToCart(currentKep, mu)';
        end
    end
end