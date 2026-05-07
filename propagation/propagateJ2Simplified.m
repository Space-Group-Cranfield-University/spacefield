function kep_next = propagateJ2Simplified(dt, kep_prev, CONST)
    if nargin < 3
        CONST = initializeAstronomicalConstants;
    end
    kep_next = kep_prev;
    T = computeOrbitalPeriod(kep_prev(1), CONST.MU_E);
    n = getMeanMotion(kep_prev(1), CONST.MU_E);
    M_prev = convert_anomaly_v2M(kep_prev(end), kep_prev(2));
    M_next = M_prev + n * mod(dt, T);
    kep_next(end) = convert_anomaly_M2v(M_next, kep_prev(2));
    kep_next(4) = mod( kep_next(4) + getRaanRateJ2(kep_prev, CONST) * dt, 2*pi);
    kep_next(5) = mod( kep_next(5) + getAopRateJ2(kep_prev, CONST) * dt, 2*pi);
end