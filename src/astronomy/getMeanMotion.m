function n = getMeanMotion(a, mu)
    if nargin < 2
        mu = initializeAstronomicalConstants().MU_E;
    end
    if max(size(a)) == 6
        n = sqrt(mu / a(1)^3);
    else
        n = sqrt(mu/a^3);
    end
end