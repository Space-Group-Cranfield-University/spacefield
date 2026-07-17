function TRG = extractTracksFromVisibilityMatrix(TRG, targetVisibilityMat)
    for k = 1:size(TRG, 2)
        tPrev = 0;
        trackIndex = 1;
        TRG(k).track = [];
        for j = 1:size(targetVisibilityMat, 2)
            if targetVisibilityMat(k, j) > 0 && tPrev == 0
                TRG(k).track(trackIndex, :) = zeros(1, 2);
                TRG(k).track(trackIndex, 1) = j;
                tPrev = 1;
            end
            if targetVisibilityMat(k, j) == 0 && tPrev == 1
                TRG(k).track(trackIndex, 2) = j;
                tPrev = 0;
                trackIndex = trackIndex + 1;
            end
        end
        if ~isempty(TRG(k).track) && TRG(k).track(end, 2) == 0
            TRG(k).track(trackIndex, 2) = size(targetVisibilityMat, 2);
        end
    end
end