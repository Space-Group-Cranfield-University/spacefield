function rMat = getConstellationPositionMatrix(Constellation, timestep)
    rMat = zeros(size(Constellation, 2), 3);
    if nargin < 2
        for k = 1:size(Constellation, 2)
            rMat(k, :) = Constellation(k).x0(1:3)';
        end
        return
    end
    for k = 1:size(Constellation, 2)
        rMat(k, :) = Constellation(k).xMat(timestep, 1:3);
    end
end