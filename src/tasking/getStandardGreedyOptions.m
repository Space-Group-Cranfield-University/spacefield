function OPTIONS = getStandardGreedyOptions(unconstrained)
    if nargin < 1
        unconstrained = 0;
    end
    OPTIONS.nSamples = 20;
    if unconstrained
        OPTIONS.constrained = 0;
        OPTIONS.maxSlewAngle = deg2rad(30);
        OPTIONS.isTargetStateUncertain = 0;
        OPTIONS.lowPower = 0;
    else
        OPTIONS.constrained = 1;
        OPTIONS.maxSlewAngle = deg2rad(5);
        OPTIONS.minEl = -deg2rad(70);
        OPTIONS.maxEl = deg2rad(5);
        OPTIONS.isTargetStateUncertain = 1;
        % If we want to idle during Eclipse:
        OPTIONS.lowPower = 1;
    end

    % Observers may retain inaccurate knowledge of target states if they
    % have not communicated with ground for long.
    OPTIONS.positionError = 10; % km
    OPTIONS.velocityError = OPTIONS.positionError * 1e-3; % km / s

    % A target that had been observed at the previous time step can be
    % valued more. This is to guarantee minimum track window and avoid fast
    % sensor switching.
    OPTIONS.valueOfPreviouslyObserved = 1;

    % Greedy may optimize across the complete network or each sensor on its
    % own. The second strategy is suboptimal but compatible with on-board
    % search of best pointing direction.
    OPTIONS.optimizeForCompleteNetwork = 0;
end