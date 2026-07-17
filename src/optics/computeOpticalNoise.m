function N = computeOpticalNoise(n_p, S, tau, CAMERA)
    if nargin < 4
        CAMERA = getStandardCamera;
    end
    K = CAMERA.K;
    I = CAMERA.I_dark;
    R_N_sq = CAMERA.R_N_sq;
    N = sqrt( S + n_p .* ( ( K + I ) .* tau + R_N_sq ) );
end