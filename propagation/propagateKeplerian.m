function xNext = propagateKeplerian(xPrev, deltaT, mu)
    if nargin < 3
        mu = initializeAstronomicalConstants().MU_E;
    end
    kepPrev = convertCartToKep(xPrev, mu);
    n = getMeanMotion(kepPrev(1), mu);
    M_prev = convert_anomaly_v2M(kepPrev(6), kepPrev(2));
    M_next = M_prev + n * deltaT;
    kepNext = kepPrev;
    kepNext(6) = convert_anomaly_M2v(M_next, kepNext(2));
    xNext = convertKepToCart(kepNext, mu);
end