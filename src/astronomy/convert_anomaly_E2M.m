% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                           function M = convert_anomaly_E2M(E, e)
%                                       Ver. 1.0
%
%           This function converts eccentric anomaly to mean anomaly.
%   
%   inputs:
%           E       eccentric anomaly [rad]
%           e       eccentricity [-]
%
%   outputs:
%           M       mean anomaly [rad]
%

function M = convert_anomaly_E2M(E, e)
    M = E - e*sin(E);
end