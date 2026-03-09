function TRG = initializeSingleTarget(kep, D_t)
    MU = initializeAstronomicalConstants().MU_E;
    TRG.a = kep(1);
    TRG.e = kep(2);
    TRG.in = kep(3);
    TRG.raan = kep(4);
    TRG.aop = kep(5);
    TRG.v = kep(6);
    TRG.kep = kep;
    TRG.x0 = convertKepToCart(TRG.kep, MU);
    TRG.D_t = D_t;
end