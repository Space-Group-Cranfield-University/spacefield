function sample = sampleFromBounds(bounds, isInteger)
    if nargin < 2
        isInteger = 0;
    end
    if isscalar(bounds)
        sample = bounds;
    elseif isInteger
        sample = bounds(1) - 1 + randi(bounds(2) - bounds(1) + 1);
    else
        sample = bounds(1) + rand*(bounds(2) - bounds(1));
    end
end