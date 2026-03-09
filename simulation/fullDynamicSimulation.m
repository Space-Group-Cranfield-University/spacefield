function [RESULTS, TRG, OBS] = ...
        fullDynamicSimulation(timeVec, OBS, TRG, dirSunMat, OPTIONS)
    % Requires proper Observer, Target, and Sun direction initialization
    % and propagation.
    if nargin < 5
        OPTIONS = getStandardSimulationOptions;
    end
    if nargin < 4
        dirSunMat = [];
    end
    if ~isfield(OPTIONS, "verbose")
        OPTIONS.verbose = getStandardSimulationOptions().verbose;
    end
    if ~isfield(OPTIONS, "taskingStrategy")
        OPTIONS.taskingStrategy = ...
            getStandardSimulationOptions().taskingStrategy;
    end
    if ~isfield(OPTIONS, "filterFlag")
        OPTIONS.filterFlag = getStandardSimulationOptions().filterFlag;
    end
    if ~isfield(OPTIONS, "postProcessing")
        OPTIONS.postProcessing = ...
            getStandardSimulationOptions().postProcessing;
    end
    if ~isfield(OPTIONS, "getCoverage")
        OPTIONS.getCoverage = getStandardSimulationOptions().getCoverage;
    end
    if ~isfield(OPTIONS, "maxOrbitalFraction")
        OPTIONS.maxOrbitalFraction = ...
            getStandardSimulationOptions().maxOrbitalFraction;
    end
    if ~isfield(OPTIONS, "minTrackLength")
        OPTIONS.minTrackLength = ...
            getStandardSimulationOptions().minTrackLength;
    end
    if ~isfield(OPTIONS, "TASKING_OPTIONS")
        OPTIONS.TASKING_OPTIONS = getStandardGreedyOptions;
    end

    deltaT = timeVec(2) - timeVec(1);
    tic
    sensorActivationMat = zeros(size(OBS, 2), size(timeVec, 2));
    targetVisibilityMat = zeros(size(TRG, 2), size(timeVec, 2));
    ObsTrgCrossVisibilityMat = zeros(size(OBS, 2), size(TRG, 2));

    for j = 2:size(timeVec, 2)
        % Task sensor (this is the most computationally expensive step if
        % using a complex tasking strategy)
        OBS = taskSensors...
            (j, OBS, OPTIONS.taskingStrategy, TRG, ...
            ObsTrgCrossVisibilityMat, dirSunMat(j, :)', ...
            OPTIONS.TASKING_OPTIONS);

        % Compute which targets are visible to which observers.
        ObsTrgCrossVisibilityMat    = fullObsTrgVisibilityTest...
                                    (j, OBS, TRG, dirSunMat);
        sensorActivationMat(:, j)   = sum(ObsTrgCrossVisibilityMat, 2);
        targetVisibilityMat(:, j)   = sum(ObsTrgCrossVisibilityMat', 2);

        % Apply Kalman filtering using the available measurements.
        % Computationally expensive, therfore it is possible to disable
        % if one does not want the estimated states and covariances but
        % only metrics such as number of tracks, track length, and so
        % on.
        if OPTIONS.filterFlag
            TRG = updateAllTargetStates...
                (j, deltaT, TRG, OBS, ObsTrgCrossVisibilityMat);
        end

        % Display step
        if OPTIONS.verbose && ~mod(j, 100)
            disp(string(j) + " / "+string(size(timeVec, 2)))
        end
    end
    RESULTS.sensorActivationMat = sensorActivationMat;
    RESULTS.targetVisibilityMat = targetVisibilityMat;
    if OPTIONS.verbose
        disp("Simulation complete, time [s]: "+toc)
    end
    tic
    if OPTIONS.postProcessing
        [TRG, ~, ~, ~, nTrackVec, ~, ~, allTrackLengths, allTbts] = ...
            computeTrackQualityMetrics(TRG, targetVisibilityMat, deltaT);
        [nTrackedTargets, nTrackedTargetsVec] = ...
            computeNumberOfTargetsTrackedOnce(targetVisibilityMat);
        [TRG, nWellTracked, nWellTrackedVec, isWellTrackedMat] = ...
            computeNumberOfWellTrackedTargets...
            (TRG, targetVisibilityMat, deltaT, ...
            OPTIONS.maxOrbitalFraction, OPTIONS.minTrackLength);
        completeness = computeCompleteness(targetVisibilityMat);
        RESULTS.nTrackVec = nTrackVec;
        RESULTS.allTrackLengths = allTrackLengths;
        RESULTS.allTbts = allTbts;
        RESULTS.nTrackedTargets = nTrackedTargets;
        RESULTS.nTrackedTargetsVec = nTrackedTargetsVec;
        RESULTS.nWellTracked = nWellTracked;
        RESULTS.nWellTrackedVec = nWellTrackedVec;
        RESULTS.isWellTrackedMat = isWellTrackedMat;
        RESULTS.completeness = completeness;
        if OPTIONS.getCoverage
            [~, ~, dailyCoverageVec] = ...
                getDailyCoverage(timeVecCoverage, OBS, dirSun0);
            RESULTS.dailyCoverageVec = dailyCoverageVec;
        end
        if OPTIONS.filterFlag
            [~, estimationErrorMat, ~, normalizedEstimationErrorVec] = ...
                computeEstimationError(timeVec, TRG);
            [meanEstimationErrorUnderTracking, ...
                meanEstimationErrorUnderTrackingVec, ...
                medianEstimationErrorUnderTracking, ...
                medianEstimationErrorUnderTrackingVec] = ...
                computeEstimationErrorUnderTracking...
                (estimationErrorMat, targetVisibilityMat);
            RESULTS.estimationErrorMat = estimationErrorMat;
            RESULTS.normalizedEstimationErrorVec = ...
                normalizedEstimationErrorVec;
            RESULTS.meanEstimationErrorUnderTrackingVec = ...
                meanEstimationErrorUnderTrackingVec;
            RESULTS.meanEstimationErrorUnderTracking = ...
                meanEstimationErrorUnderTracking;
            RESULTS.medianEstimationErrorUnderTrackingVec = ...
                medianEstimationErrorUnderTrackingVec;
            RESULTS.medianEstimationErrorUnderTracking = ...
                medianEstimationErrorUnderTracking;
        end
    end
    if OPTIONS.verbose
        disp("Post-processing complete, time [s]: "+toc)
    end
end