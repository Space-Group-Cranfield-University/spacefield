% convertKepToCart converts Keplerian parameters to cartesian coordinates.
%
%   cart = convertKepToCart(kep, mu)
%
%   Inputs:
%   - kep   : Keplerian parameters. [a, e, in, raan, aop, v]
%   - mu    : gravitational parameter of the central body. [km^3/s^2]
%
%   Outputs:
%   - cart  : 6x1 vector containing cartesian coordinates (position and
%           velocity).

function cart = convertKepToCart(kep, mu)

    a = kep(1); % semi-major axis [km]
    e = kep(2); % eccentricity [-]
    i = kep(3); % inclination [rad]
    OM = kep(4); % RAAN [rad]
    om = kep(5); % aop [rad]
    theta = kep(6); % True anomaly [rad]

    % Compute eccentric anomaly
    E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));

    % Compute distance from central body
    rnorm = a*(1-e*cos(E));

    % Compute velocity coefficient
    vcoeff = sqrt(mu*a)/rnorm;

    % Compute r, v in orbital frame
    r_o = rnorm*[cos(theta) sin(theta) 0]';
    v_o = vcoeff*[-sin(E) sqrt(1-e^2)*cos(E) 0]';

    % Rotation matrices
    R_x = @(t) [1 0 0;
                0 cos(t) sin(t);
                0 -sin(t) cos(t)];
    R_z = @(t) [cos(t) sin(t) 0;
                -sin(t) cos(t) 0;
                0 0 1];

    % Compute rotation matrix
    R = R_z(-OM)*R_x(-i)*R_z(-om);

    % Compute cartesian coordinates
    r = R*r_o;
    v = R*v_o;
    cart = [r;v];

end