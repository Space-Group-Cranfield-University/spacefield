function dom = getAopRateJ2(kep, CONST)
    % https://conference.sdo.esoc.esa.int/proceedings/sdc7/paper/530/SDC7-paper530.pdf
    if nargin < 2
        CONST =  initializeAstronomicalConstants;
    end
    n = getMeanMotion(kep(1), CONST.MU_E);
    dom = 3 * (1 - 5 / 4 * sin(kep(3))^2 ) * CONST.J2 * CONST.R_E_NORAD^2 * n / ( kep(1) * (1 - kep(2)^2) )^2;
end