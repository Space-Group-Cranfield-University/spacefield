% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                           function E = convert_anomaly_v2E(v, e)
%                                       Ver. 1.0
%
%           This function converts true anomaly to eccentric anomaly.
%   
%   inputs:
%           v       true anomaly [rad]
%           e       eccentricity [-]
%
%   outputs:
%           E       eccentric anomaly [rad]
%

function E = convert_anomaly_v2E(v, e)
    E = 2*atan(sqrt((1-e)/(1+e))*tan(v/2));
end