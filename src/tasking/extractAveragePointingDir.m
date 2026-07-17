function [meanAz, meanEl, stdAz, stdEl] = extractAveragePointingDir(timeVec, OBS)
    % Mean
    meanAz = zeros(1, size(timeVec, 2));
    meanEl = zeros(1, size(timeVec, 2));
    for k = 1:size(OBS, 2)
        for j = 1:size(timeVec, 2)
            meanAz(j) = meanAz(j) + OBS(k).pointingAzElMat(j, 1);
            meanEl(j) = meanEl(j) + OBS(k).pointingAzElMat(j, 2);
        end
    end 
    meanAz = meanAz / size(OBS, 2);
    meanEl = meanEl / size(OBS, 2);
    % Standard deviation
    varAz = zeros(1, size(timeVec, 2));
    varEl = zeros(1, size(timeVec, 2));
    for k = 1:size(OBS, 2)
        for j = 1:size(timeVec, 2)
            varAz(j) = varAz(j) + ...
                ( OBS(k).pointingAzElMat(j, 1) - meanAz(j) )^2;
            varEl(j) = varEl(j) + ...
                ( OBS(k).pointingAzElMat(j, 2) - meanEl(j) )^2;
        end
    end
    varAz = varAz / ( size(OBS, 2) - 1 );
    varEl = varEl / ( size(OBS, 2) - 1 );
    varTh = varAz + varEl;
    stdAz = sqrt( varAz );
    stdEl = sqrt( varEl );
end