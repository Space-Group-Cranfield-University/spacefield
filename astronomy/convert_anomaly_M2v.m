% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                           function v = convert_anomaly_M2v(M, e)
%                                       Ver. 1.0
%
%           This function converts mean anomaly to true anomaly.
%   
%   inputs:
%           M       mean anomaly [rad]
%           e       eccentricity [-]
%
%   outputs:
%           v       true anomaly [rad] 
%

function v = convert_anomaly_M2v(M, e)
    v = convert_anomaly_E2v(convert_anomaly_M2E(M, e), e);
end