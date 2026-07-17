function currentObservations = getCurrentObservations(ObsTrgCrossVisibilityMat)
    % Converts the observer-target cross-visibility matrix (at a single
    % time step) to a list of observations (for that time step) of the form
    % currentObservation(k) = (OBS_ID, TRG_ID), k is the observation number
    nonZeroIndexesObs = find(ObsTrgCrossVisibilityMat == 1);
    nObs = size(ObsTrgCrossVisibilityMat, 1);
    obsID = mod(nonZeroIndexesObs, nObs);
    obsID(obsID == 0) = nObs;
    trgID = 1 + (nonZeroIndexesObs - obsID) / nObs;
    currentObservations = [obsID, trgID];
end