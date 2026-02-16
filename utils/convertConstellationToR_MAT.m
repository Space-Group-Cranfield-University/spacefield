function rMat = convertConstellationToR_MAT(Constellation)
    rMat = [];
    for k = 1:size(Constellation, 2)
        rMat = [rMat; Constellation(k).x0(1:3)'];
    end
end