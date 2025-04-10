% initialiseAstronomicalConstants Initialises useful constants.
%
%   Constants = initialiseAstronomicalConstants() returns a structure
%   Constants containing the following fields (alphabetical order):
%   .C              : Speed of light [m/s]
%   .ECL            : Earth equator-ecliptic inclination [rad]
%   .G_0            : Standard Earth gravity [m/s^2]
%   .H_ATM          : Atmosphere scale height [m]
%   .H_PLANCK       : Planck's constant [J*s]
%   .H_SOLAR        : Solar constant [W/m^2]
%   .J2             : Earth spherical harmonic degree 2 [-]
%   .LAMBDA_VISIBLE : Average wavelength visible light [m]
%   .M_S            : Apparent magnitude of the Sun [-]
%   .MU_E           : Earth gravitational parameter [km^3/s^2]
%   .NI_VISIBLE     : Average frequency visible light [Hz]
%   .OM_E           : Earth spin rate [rad/s]
%   .OM_S           : Earth revolution rate around the Sun [rad/s]
%   .R_E            : Average Earth radius [km]
%   .R_E_NORAD      : NORAD Earth radius [km]
%   .R_S            : Average Sun-Earth distance / 1 Astronomical Unit [km]
%   .RHO_NORAD      : NORAD atmospheric density [kg/(m^2*R_E_NORAD)]
%   .RHO_0          : ISA air density at standard conditions [kg/m^3]

function Constants = initialiseAstronomicalConstants()

    % Units are km, kg, s, rad. Bstar is 1/Re and has std value 0.0002
    Constants.C = 299792458; % speed of light [m/s]
    Constants.ECL = deg2rad(23.4); % ecliptic angle [rad]
    Constants.G_0 = 9.80665; % standard Earth's gravity [m/s^2]
    Constants.H_ATM = 8.5*1e3; % atmosphere scale height [m]
    Constants.H_PLANCK = 6.62607015*1e-34; % Planck's constant [J*s]
    Constants.H_SOLAR = 1361; % Solar constant [W/m^2]
    Constants.J2 = 1.082e-3; % Earth's spherical harmonic degree 2 [-]

    % average wavelength visible light [m]
    Constants.LAMBDA_VISIBLE = 550*1e-9;

    % Earth's gravitational parameter [km^3/s^2]
    Constants.MU_E = 398600.5;

    % average frequency visible light [Hz]
    Constants.NI_VISIBLE = Constants.C / Constants.LAMBDA_VISIBLE;

    % Earth's rotation rate [rad/s]
    Constants.OM_E = 2*pi/(3600*23+60*56);

    % Earth's revolution rate around Sun [rad/s]
    Constants.OM_S = 2*pi/(3600*24*365+3600*6+60*9);

    Constants.R_E = 6371; % Average Earth radius [km]
    Constants.R_E_NORAD = 6378.135; % NORAD Earth radius [km]

    % Sun-Earth distance approximated with Earth heliocentric semi-major 
    % axis (circular orbit assumption) [km]
    Constants.R_S = 149598023;

    % NORAD atmospheric density [kg/(m^2*Re)]
    Constants.RHO_NORAD = 2.461e-5*6378.135;

    % International Standard Atmosphere dry air density at standard 
    % conditions [kg/m^3]
    Constants.RHO_0 = 1.225;

    Constants.M_S = -26.73;

end