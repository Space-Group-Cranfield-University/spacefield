% ------------------------------------------------------------------------------
%           Author: Antonio D'Anniballe
%
%                           function E = convert_anomaly_M2E(M, e)
%                                       Ver. 0.1
%
%           This function converts mean anomaly to eccentric anomaly.
%   
%   inputs:
%           M       mean anomaly [rad]
%           e       eccentricity [-]
%
%   outputs:
%           E       eccentric anomaly [rad]
%

function E = convert_anomaly_M2E(M, e)
    func = @(E) E - e*sin(E) - M;
    E = fsolve(func, M, optimset('Display','off'));
end