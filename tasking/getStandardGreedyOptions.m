function OPTIONS = getStandardGreedyOptions()
    OPTIONS.constrained = 1;
    if OPTIONS.constrained
        OPTIONS.maxSlewAngle = deg2rad(1);
    else
        OPTIONS.maxSlewAngle = deg2rad(30);
    end
    OPTIONS.minEl = deg2rad(0);
    OPTIONS.maxEl = deg2rad(10);
    OPTIONS.nSamples = 100;

    % Observers may retain inaccurate knowledge of target states if they
    % have not communicated with ground for long.
    OPTIONS.positionError = 0; % km
    OPTIONS.velocityError = 0; % km / s

    % A target that had been observed at the previous time step can be
    % valued more. This is to guarantee minimum track window and avoid fast
    % sensor switching.
    OPTIONS.valueOfPreviouslyObserved = 1;

    % Greedy may optimize across the complete network or each sensor on its
    % own. The second strategy is suboptimal but compatible with on-board
    % search of best pointing direction.
    OPTIONS.optimizeForCompleteNetwork = 1;
end