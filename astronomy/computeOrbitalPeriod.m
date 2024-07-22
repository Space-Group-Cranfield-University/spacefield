% computeOrbitalPeriod computes the period of a given orbit around a given
% central body.
%
%   orbitalPeriod = computeOrbitalPeriod(a, mu)
%
%   Inputs:
%   - a     : orbital semi-major axis [km].
%   - mu    : gravitational parameter of the central body [km^3/s^2].
%
%   Outputs:
%   - orbitalPeriod : orbital period [s].

function orbitalPeriod = computeOrbitalPeriod(a, mu)

    orbitalPeriod = 2*pi / sqrt(mu) * a^1.5;

end