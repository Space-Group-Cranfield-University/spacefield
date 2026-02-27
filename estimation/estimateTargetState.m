function TRG = estimateTargetState(deltaT, timestep, kTarget, TRG, y, h)
    f = @(x) propagateKeplerian(x, deltaT);
    xPrev = TRG(kTarget).xMatEst(timestep - 1, :)';
    P_prev = TRG(kTarget).P(:, :, timestep - 1);
    if ~isempty(y)
        [xNext, P_next] = stepUKF(xPrev, P_prev, y, f, h);
    else
        [xNext, P_next] = predictUT(xPrev, P_prev, f);
    end
    TRG(kTarget).xMatEst(timestep, :) = xNext';
    TRG(kTarget).P(:, :, timestep) = P_next;
end