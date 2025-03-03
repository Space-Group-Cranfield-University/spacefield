% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                               function kep = convert_cart2kep(y, mu)
%                                       Ver. 0.1
%
%           This function transforms a cartesian state vector into its
%           corresponding keplerian elements.
%   
%   inputs:
%           y       cartesian state vector as column vector [km] and [km/s]
%           mu      gravitational parameter [km^3/s^2]
%
%   outputs:
%           kep     keplerian state vector [km, -, rad, rad, rad, rad]
%
%   notes:
%           This function needs further testing! Careful when using it with
%           singular orbits. This function currently handles ONLY true
%           anomaly. Mean and eccentric anomaly to be added.
%

function kep = convertCart2Kep(y, mu)
    % WARNING: Function does NOT handle singularity at i = 90 deg
    kep = zeros(6,1);
    r = y(1:3,1);
    v = y(4:6,1);
    R_x = @(t) [1 0 0;
                0 cos(t) sin(t);
                0 -sin(t) cos(t)]; % rotation matrix around local x
    R_z = @(t) [cos(t) sin(t) 0;
                -sin(t) cos(t) 0;
                0 0 1]; % rotation matrix around local z
    if norm(r) < 1e-15
        disp("Error in convert_cart2kep() - Null radial distance")
        return
    end
    a = inv(2/norm(r)-norm(v)^2/mu);
    h = cross(r,v);
    if norm(h) < 1e-15
        disp("Error in convert_cart2kep() - Null angular momentum")
        return
    end   
    i = acos(h(3)/norm(h));
    e_vect = cross(v,h)/mu - r/norm(r);
    e = norm(e_vect);
    N = cross([0 0 1], h);
    if e < 1e-15
        disp("Warning in convert_cart2kep() - Circular orbit")
        om = 0;
        if norm(N) < 1e-15
            disp("Warning in convert_cart2kep() - Equatorial orbit")
            OM = 0;
            theta = atan2(r(2),r(1));
            kep = [a; e; i; OM; om; theta];
            return
        end
        if N(2) < 0
            OM = 2*pi - acos(N(1)/norm(N));
        else
            OM = acos(N(1)/norm(N));
        end
        r_o = R_x(i)*R_z(OM)*r; % rotation into orbital frame
        theta = atan2(r_o(2),r_o(1));
        kep = [a; e; i; OM; om; theta];
        return
    end
    if r'*v < 0
        theta = 2*pi - acos(e_vect'*r/(e*norm(r)));
    else
        theta = acos(e_vect'*r/(e*norm(r)));
    end
    if norm(N) < 1e-15
        disp("Warning in convert_cart2kep() - Equatorial orbit")
        OM = 0;
        lon = atan2(r(2),r(1));
        om = lon-theta;
        kep = [a;e;i;OM;om;theta];
        return
    end
    if e_vect(3) < 0
        om = 2*pi - acos(N*e_vect/(e*norm(N)));
    else
        om = acos(N*e_vect/(e*norm(N)));
    end
    if N(2) < 0
        OM = 2*pi - acos(N(1)/norm(N));
    else
        OM = acos(N(1)/norm(N));
    end
    kep = [a;e;i;OM;om;theta];
end