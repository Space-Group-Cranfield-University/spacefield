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
%   NOTICE: ONLY works for circular (or quasi-circular) constellation orbits

function Constellation = propagateConstellation(timeVec, Constellation, includeJ2, CONST)
    if nargin < 3
        includeJ2 = 0;
    end
    if nargin < 5
        CONST = initializeAstronomicalConstants;
    end
    if includeJ2
        Constellation = propagateConstellationJ2(timeVec, Constellation, CONST);
    else
        Constellation = propagateConstellationKeplerian(timeVec, Constellation, CONST.MU_E);
    end
end

function Constellation = propagateConstellationKeplerian(timeVec, Constellation, mu)
    for iSat = 1:size(Constellation, 2)
        T = Constellation(iSat).T;
        thetaVec = 2*pi * timeVec / T;
        Constellation(iSat).xMat = zeros(size(thetaVec, 2), 6);
        for iTheta = 1:size(thetaVec, 2) 
            currentKep = [Constellation(iSat).kep(1:5, 1); Constellation(iSat).kep(6, 1) + thetaVec(iTheta)];
            Constellation(iSat).xMat(iTheta, :) = convertKepToCart(currentKep, mu)';
        end
    end
end

function Constellation = propagateConstellationJ2(timeVec, Constellation, CONST) 
    dt = timeVec(2) - timeVec(1);
    for iSat = 1:size(Constellation, 2)
        Constellation(iSat).xMat = zeros(size(timeVec, 2), 6);
        in = Constellation(iSat).in;
        a = Constellation(iSat).a;
        T = Constellation(iSat).T;
        thetaVec = 2*pi * timeVec / T;
        n = 2 * pi / T;
        for iTheta = 1:size(thetaVec, 2) 
            currentKep = [Constellation(iSat).kep(1:5, 1); Constellation(iSat).kep(6, 1) + thetaVec(iTheta)];
            Om = currentKep(5);
            dOm = getRaanRateJ2(currentKep, CONST);
            Om = Om + dOm * dt;
            currentKep(5) = Om;
            Constellation(iSat).xMat(iTheta, :) = convertKepToCart(currentKep, CONST.MU_E)';
        end
    end
end