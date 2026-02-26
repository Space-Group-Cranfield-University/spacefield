function a = getPerturbationDrag(x, B, densityModel, CONST)
    if nargin < 4
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 3
        densityModel = "NRLMSIS-00";
    end
    if nargin < 2
        B = 0.01;
    end
    rho = 1e3 * getAtmosphericDensity(x, densityModel, CONST);
    a = - 0.5 * B * rho * norm(x(4:6)) * x(4:6, 1);
end