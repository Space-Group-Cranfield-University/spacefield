% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                           function v = convert_anomaly_E2v(E, e)
%                                       Ver. 1.0
%
%           This function converts eccentric anomaly to true anomaly.
%   
%   inputs:
%           E       eccentric anomaly [rad]
%           e       eccentricity [-]
%
%   outputs:
%           v       true anomaly [rad] 
%

function v = convert_anomaly_E2v(E, e)
    v = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end