function [nSatPareto, coveragePareto] = estimateParetoFrontier(nSatVec, coverageVec)
    [nSatVecSorted, sortIndex] = sort(nSatVec);
    coverageVecSorted = coverageVec(sortIndex);
    coveragePareto = coverageVecSorted(1);
    nSatPareto = nSatVecSorted(1);
    for k = 2:length(coverageVecSorted)
        if nSatVecSorted(k) == nSatPareto(end) % Case where two coverage values have the same number of satellites
            if coverageVecSorted(k) > coveragePareto(end)
                coveragePareto(end) = coverageVecSorted(k);
                nSatPareto(end) = nSatVecSorted(k);
            end
        elseif coverageVecSorted(k) > coveragePareto(end)
            coveragePareto = [coveragePareto, coverageVecSorted(k)];
            nSatPareto = [nSatPareto, nSatVecSorted(k)];
        end
    end
end