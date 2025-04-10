% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                           function M = convert_anomaly_v2M(v, e)
%                                       Ver. 1.0
%
%           This function converts true anomaly to mean anomaly.
%   
%   inputs:
%           v       true anomaly [rad]
%           e       eccentricity [-]
%
%   outputs:
%           M       mean anomaly [rad]
%

function M = convert_anomaly_v2M(v, e)
    M = convert_anomaly_E2M(convert_anomaly_v2E(v,e),e);
end