function V = computeApparentMagnitude(deltaR, A_t, alpha, TRG_DATA)
    % Computes apparent magnitude of a given target.
    %
    % V = computeApparentMagnitude(deltaR, A_t, alpha, TRG_DATA)
    % V = computeApparentMagnitude(deltaR, A_t, alpha)
    % V = computeApparentMagnitude(deltaR, A_t)
    if nargin < 4
        TRG_DATA = getStandardTargetReflectivityData;
    end
    if nargin < 3
        alpha = 0;
    end
    Phi = computeMixedPhaseFunction(alpha, TRG_DATA);
    V = - 26.74 - 2.5 * log10( A_t * Phi ) + 5 * log10( deltaR );
end