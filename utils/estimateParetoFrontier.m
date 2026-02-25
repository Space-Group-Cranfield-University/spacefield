function [nSatParetoAugmented, coverageParetoAugmented] = estimateParetoFrontier(nSatVec, coverageVec, costOrFitness)
    if nargin < 3
        costOrFitness = "CostFitness";
    end
    if strcmp(costOrFitness, "FitnessFitness")
        [nSatParetoAugmented, coverageParetoAugmented] = estimateParetoFrontierFitnessFitness(nSatVec, coverageVec);
    else
        [nSatParetoAugmented, coverageParetoAugmented] = estimateParetoFrontierCostFitness(nSatVec, coverageVec);
    end
end

function [nSatPareto, coveragePareto] = estimateParetoFrontierFitnessFitness(nSatVec, coverageVec)
    [nSatVecSorted, sortIndex] = sort(nSatVec, 'descend');
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

function [nSatParetoAugmented, coverageParetoAugmented] = estimateParetoFrontierCostFitness(nSatVec, coverageVec)
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
    nSatParetoAugmented = min(nSatPareto):max(nSatPareto);
    coverageParetoAugmented = zeros(size(nSatVec));
    indexPareto = 1;
    for j = 1:size(nSatParetoAugmented, 2)
        nSat = nSatParetoAugmented(j);
        if nSat >= nSatPareto(indexPareto)
            indexPareto = indexPareto + 1;
        end
        coverageParetoAugmented(j) = coveragePareto(indexPareto-1);
    end
end