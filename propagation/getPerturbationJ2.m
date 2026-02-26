function a = getPerturbationJ2(x, CONST)
    if nargin < 2
        CONST = initialiseAstronomicalConstants;
    end
    R_e = CONST.R_E_NORAD;
    rNorm = norm(x(1:3));
    X = x(1);
    y = x(2);
    z = x(3);
    k = - 1.5 * CONST.J2 * CONST.MU_E * R_e^2 / rNorm^5;
    a = k * ...
        [ X * (1 - 5*z^2/rNorm^2) ; ...
          y * (1 - 5*z^2/rNorm^2) ; ...
          z * (3 - 5*z^2/rNorm^2) ];
end