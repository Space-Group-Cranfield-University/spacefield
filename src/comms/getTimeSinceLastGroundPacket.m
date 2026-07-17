function timeSinceLastGroundPacket = getTimeSinceLastGroundPacket(timeVec, OBS, stationMat, OPTIONS, CONST)
    if nargin < 5
        CONST = initializeAstronomicalConstants;
    end
    % Requires propagated observers and station locations
    deltaT = timeVec(2) - timeVec(1);
    timeSinceLastGroundPacket = zeros(size(timeVec, 2), size(OBS, 2));
    for j = 2:size(timeVec, 2)
        DCM = getRotationMatrixECEFToECI(timeVec(j));
        % Satellite-to-ground communication
        for k = 1:size(OBS, 2)
            timeSinceLastGroundPacket(j, k) = timeSinceLastGroundPacket(j-1, k) + deltaT;
            rObs = OBS(k).xMat(j, 1:3)';
            for i = 1:size(stationMat, 1)
                uECEF = convertRaDecToPosition([1; stationMat(i, 1:2)']);
                uECI = DCM * uECEF;
                rECI = CONST.R_E * uECI;
                if acos(uECI' * (rObs - rECI) / norm(rObs - rECI)) < (pi/2 - OPTIONS.minEl)
                    timeSinceLastGroundPacket(j, k) = 0;
                end
            end
        end
        % Satellite-to-satellite communication
        if OPTIONS.satToSatCommsEnabled
            for k = 1:size(OBS, 2)
                rObs_1 = OBS(k).xMat(j, 1:3)';
                for i = 1:(size(OBS, 2)-1)
                    index = i;
                    if i >= k
                        index = index + 1;
                    end
                    currTime = timeSinceLastGroundPacket(j, k);
                    advTime = timeSinceLastGroundPacket(j, index);
                    rObs_2 = OBS(index).xMat(j, 1:3)';
                    dr = norm(rObs_2 - rObs_1);
                    if advTime < currTime && dr < OPTIONS.maxDistance && isAboveTheHorizon(rObs_1, rObs_2, 0, 6371+100)
                        timeSinceLastGroundPacket(j, k) = advTime;
                    end
                end
            end
        end
    end
end