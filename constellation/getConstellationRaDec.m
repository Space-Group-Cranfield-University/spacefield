function OBS = getConstellationRaDec(OBS)
    for j = 1:size(OBS, 2)
        OBS(j).raDecMat = zeros(size(OBS(j).xMat, 1), 3);
        for k = 1:size(OBS(j).xMat, 1)          
            r = OBS(j).xMat(k, :);
            raDec = convertPositionToRaDec(r);
            OBS(j).raDecMat(k, :) = raDec';
        end
    end
end