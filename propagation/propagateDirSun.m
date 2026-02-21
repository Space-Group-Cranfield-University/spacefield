function dirSunMat = propagateDirSun(timeVec, dirSun0, Constants)
    % Simple function for propagating Sun direction (circular orbit assumption)
    % ECL   : ecliptic frame (inertial, ecliptic plane)
    % HCI   : heliocentric inertial frame (inertial, equatorial plane)
    % ECI   : Earth-centred inertial frame
    if nargin < 3
        Constants = initializeAstronomicalConstants;
    end

    th = - Constants.ECL; % Ecliptic angle (Earth-Sun is -23.4°)
    N_t = size(timeVec, 2);

    % DCM rotating ECL to HCI
    dcmEclToHci = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];

    dirSunMat = nan(N_t, 3);
    dirEarthEcl0 = - dcmEclToHci' * dirSun0;
    if abs(dirEarthEcl0(3)) > 1e-5
        warning("Sun does not lie on ecliptic.")
    end
    theta0 = atan2(dirEarthEcl0(2),dirEarthEcl0(1));
    for k = 1:N_t
        % Earth direction in ECL
        dirEarthEclX = cos(theta0 + Constants.OM_S*timeVec(k));
        dirEarthEclY = sin(theta0 + Constants.OM_S*timeVec(k));
        dirEarthEcl = [dirEarthEclX dirEarthEclY 0]';

        % Sun direction in ECI
         dirSunTemp = - dcmEclToHci * dirEarthEcl;
         dirSunMat(k,:) = dirSunTemp';
    end
end