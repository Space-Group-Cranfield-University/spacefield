function Q = getProcessNoiseCovariance(deltaT, sigmaAcceleration)
    if nargin < 2
        sigmaAcceleration = 1e-6; % compatible with high LEO drag
    end
    Q = sigmaAcceleration^2 * ...
        [deltaT^4*eye(3)/4, deltaT^3*eye(3)/2; ...
        deltaT^3*eye(3)/2, deltaT^2*eye(3)];
end