function OBS = setInitialPointingDirections(timeVec, OBS, raDec0, pointingFrame)
    % Here we assume z in the inward radial direction and x in the
    % positive out-of-plane direction for the Local frame. Other options
    % may be LVLH and ECI.
    if nargin < 4
        pointingFrame = "Local";
    end
    if nargin < 3
        raDec0 = [0, 0]';
    end
    pointingDir = convertRaDecToPosition(1, raDec0(1), raDec0(2));
    for k = 1:size(OBS, 2)
        OBS(k).pointingMat = zeros(size(timeVec, 2), 3);
        OBS(k).pointingMat(:, 1) = pointingDir(1);
        OBS(k).pointingMat(:, 2) = pointingDir(2);
        OBS(k).pointingMat(:, 3) = pointingDir(3);
        OBS(k).pointingFrame = pointingFrame;
    end
end