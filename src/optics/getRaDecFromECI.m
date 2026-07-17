function raDec = getRaDecFromECI(x_ECI)
    % Results in radians
    rNorm = norm(x_ECI(1:3));
    dec = asin(x_ECI(3)/rNorm);
    ra = (1-sign(x_ECI(2)))*pi + atan2(x_ECI(2), x_ECI(1));
    raDec = [ra; dec];
end