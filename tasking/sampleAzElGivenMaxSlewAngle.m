function [deltaAzElMat, newAzElMat] = sampleAzElGivenMaxSlewAngle(prevAzEl, OPTIONS, nSamples)
    prevEl = prevAzEl(2);
    a = 1/sqrt(cos(prevEl));
    if nargin < 3
        nSamples = OPTIONS.nSamples;
    end
    radiusSamples = OPTIONS.maxSlewAngle * sqrt(rand([1, nSamples]));
    thetaSamples = 2 * pi * rand([1,nSamples]);
    xSamples = radiusSamples .* cos(thetaSamples);
    deltaElSamples = radiusSamples .* sin(thetaSamples);
    deltaAzSamples = a * xSamples;
    index_high = find(prevEl + deltaElSamples > OPTIONS.maxEl);
    index_low = find(prevEl + deltaElSamples < OPTIONS.minEl);
    deltaElSamples(index_high) = - OPTIONS.maxSlewAngle / 2;
    deltaElSamples(index_low) = OPTIONS.maxSlewAngle / 2;
    deltaAzElMat = [deltaAzSamples; deltaElSamples];
    newAzElMat = prevAzEl + deltaAzElMat;
end