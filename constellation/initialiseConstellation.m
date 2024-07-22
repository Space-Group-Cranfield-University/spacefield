% initialiseConstellation initialises a Walker-type constellation.
%
%   Constellation = initialiseConstellation(Parameters, Constants)
%
%   Inputs:
%   - Parameters    : structure containing the following fields:
%     .nOrb         : number of constellation orbits.
%     .nSatOrb      : number of satellites per orbit.
%     .h            : constellation altitude.
%     .in           : constellation inclination.
%     .raan         : right ascension of first constellation orbit.
%     .dtheta       : phase angle between successive constellation orbits.
%   - Constants     : structure containing astronomical constants.
%
%   Outputs:
%   - Constellation : vector of structures containing satellite
%                   information. The initial state of satellite k of the
%                   constellation can be found in Constellation(k).x0.

function Constellation = initialiseConstellation(Parameters, Constants)
    
    % Extract parameters
    nOrb = Parameters.nOrb;
    nSatOrb = Parameters.nSatOrb;
    h = Parameters.h;
    in = Parameters.in;
    firstRaan = Parameters.raan;
    dtheta = Parameters.dtheta;
    
    % Initialise structure
    N_t = nSatOrb*nOrb;
    Constellation = struct('a', cell(1, N_t), 'e', cell(1, N_t),'in', cell(1, N_t),...
        'raan',cell(1, N_t),'aop', cell(1, N_t),'v', cell(1, N_t),...
        'kep',cell(1, N_t),'x0', cell(1, N_t));

    % Initialise Constellation
    for j = 1:nOrb
        raan = firstRaan + 2*pi/nOrb*(j-1);
        N = (j-1)*nSatOrb;
        for k = 1:nSatOrb
            Constellation(k+N).a = Constants.R_E + h;
            Constellation(k+N).e = 1e-6;
            Constellation(k+N).in = in;
            Constellation(k+N).raan = raan;
            Constellation(k+N).aop = 0;
            Constellation(k+N).v = (k - 1) * 2*pi/nSatOrb + (j-1) * dtheta;
            Constellation(k+N).kep  = [Constellation(k+N).a Constellation(k+N).e ...
                                    Constellation(k+N).in Constellation(k+N).raan ...
                                    Constellation(k+N).aop Constellation(k+N).v]';
            Constellation(k+N).x0 = convertKepToCart(Constellation(k+N).kep, Constants.MU_E);
        end
    end
end