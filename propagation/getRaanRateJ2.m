function dOm = getRaanRateJ2(kep, CONST)
    if nargin < 2
        CONST =  initializeAstronomicalConstants;
    end
    n = getMeanMotion(kep(1), CONST.MU_E);
    dOm = - 1.5 * CONST.J2 * CONST.R_E_NORAD^2 * n * cos(kep(3)) / kep(1)^2;
end