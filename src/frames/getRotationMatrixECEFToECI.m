function DCM = getRotationMatrixECEFToECI(time, CONST)
    if nargin < 2
        CONST = initializeAstronomicalConstants;
    end
    OM = CONST.OM_E;
    theta = OM * time;
    while theta >= 2*pi
        theta = theta - 2*pi;
    end
    DCM = [ cos(theta), -sin(theta),0; 
            sin(theta), cos(theta), 0; 
            0,          0,          1 ];
end