function a = getKeplerianAcceleration(x, CONST)
    if nargin < 2
        CONST = initializeAstronomicalConstants;
    end
    a = - CONST.MU_E / norm(x(1:3))^3 .* x(1:3, 1);
end