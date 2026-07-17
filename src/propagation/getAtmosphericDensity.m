function rho = getAtmosphericDensity(x, densityModel, CONST)
    warning('off', 'all') % Necessary to speed up density computation
    % May require Aerospace Toolbox
    if nargin < 2
        densityModel = "NRLMSIS-00";
    end
    if nargin < 3
        CONST = initializeAstronomicalConstants;
    end
    switch densityModel
        case "ISA"
            h = norm(x(1:3)) - CONST.R_E;
            rho = 101325 * exp(-0.00012 * h) / (287.05 * (273.15 - 0.0065 * h));
        case "USStandard"
            h = norm(x(1:3)) - CONST.R_E;
            rho = 1.225 * (1 - 0.0000225577 * h)^4.2561;
        case "NRLMSIS-00"
            rNorm = norm(x(1:3));
            h = rNorm - CONST.R_E;
            lat = asin(x(3) / rNorm);
            if h < 800
                [~, rho] = atmosnrlmsise00(h * 1e3, rad2deg(lat), 0, 2025, 1, 0);
                rho = rho(6);
            else
                rho = 1e-18;
            end
        otherwise
            h = norm(x(1:3)) - CONST.R_E;
            rho = CONST.RHO_0 * exp(- h / CONST.H_ATM * 1e3);
    end
end