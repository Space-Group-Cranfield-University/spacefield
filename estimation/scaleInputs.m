function [xScaled, P_scaled] = scaleInputs(x, P, D)
    if nargin < 3
        D = getStandardScaleMatrix;
    end
    xScaled = D * x;
    P_scaled = D * P * D';
end