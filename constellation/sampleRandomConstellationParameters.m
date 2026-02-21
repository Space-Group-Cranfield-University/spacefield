function Parameters = ...
            sampleRandomConstellationParameters...
            (nSatBounds, inBounds, hBounds, excludeSingleOrbits)
    isPrime = 1;
    while isPrime
        nSat = sampleFromBounds(nSatBounds, 1);
        isPrime = isprime(nSat);
    end
    in = sampleFromBounds(inBounds);
    h = sampleFromBounds(hBounds);
    if excludeSingleOrbits
        factorsOfNsat = factor(nSat);
    else
        factorsOfNsat = [ 1, factor(nSat)];
    end
    factorSubsetSize = randi(size(factorsOfNsat, 2));
    nOrbFactors = randsample(factorsOfNsat, factorSubsetSize);
    nOrb = prod(nOrbFactors);
    nSatOrb = nSat / nOrb;
    Parameters.h = h;
    Parameters.nOrb = nOrb;
    Parameters.nSatOrb = nSatOrb;
    Parameters.in = in;
    Parameters.F = randi(nOrb)-1;
end