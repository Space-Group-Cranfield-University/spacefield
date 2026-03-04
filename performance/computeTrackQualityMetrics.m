function [TRG, nTrack, trackLength, tbt, nTrackVec, trackLengthVec, tbtVec, ...
    allTrackLengths, allTbts] = ...
            computeTrackQualityMetrics(TRG, targetVisibilityMat, deltaT)
    % tbt is time between tracks
    TRG = extractTracksFromVisibilityMatrix(TRG, targetVisibilityMat);
    allTrackLengths = [];
    allTbts = [];
    for k = 1:size(TRG, 2)
        nTrackVec(k) = size(TRG(k).track, 1);
        currentTrackLength = 0;
        current_tbt = 0;
        for j = 1:nTrackVec(k)
            currentTrackLength(j) = ...
                TRG(k).track(j, 2) - TRG(k).track(j, 1);
            allTrackLengths = [allTrackLengths; currentTrackLength(j)];
            if j ~= nTrackVec(k)
                current_tbt(j) = TRG(k).track(j+1, 1) - TRG(k).track(j, 2);
                allTbts = [allTbts; current_tbt(j)];
            end
        end
        trackLengthVec(k) = mean(currentTrackLength);
        tbtVec(k) = mean(current_tbt);
    end
    nTrack = floor(mean(nTrackVec));
    trackLength = mean(trackLengthVec) * deltaT;
    tbt = mean(tbtVec) * deltaT;
end