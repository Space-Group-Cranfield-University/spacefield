function SNR = computeOpticalSNR(tau, deltaR, A_t, v_t, alpha, CAMERA, TRG_DATA, CONST)

    % Function to get optical signal-to-noise ratio for a target tracing a
    % streak over an observing detector.
    %
    % REFERENCE:
    %   - Pineau D & Felicetti L (2023) Design of an optical system for a 
    %   Multi-CubeSats debris surveillance mission, Acta Astronautica, 
    %   210 (September) 535-546.
    %
    % Assumption: measurements are taken over the interval [t_0, t_0 + tau]
    %
    % Notice: to get the result in decibel, use 10*log10(SNR).
    %
    % SNR = getOpticalSNR(tau, deltaR, A_t, v_t, alpha, CAMERA, TRG_DATA, CONST)
    % SNR = getOpticalSNR(tau, deltaR, A_t, v_t, alpha, CAMERA, TRG_DATA)
    % SNR = getOpticalSNR(tau, deltaR, A_t, v_t, alpha, CAMERA)
    % SNR = getOpticalSNR(tau, deltaR, A_t, v_t, alpha)
    % SNR = getOpticalSNR(tau, deltaR, A_t, v_t)
    % SNR = getOpticalSNR(tau, deltaR, A_t)

    if nargin < 8
        CONST = initializeAstronomicalConstants;
    end

    if nargin < 7
        TRG_DATA = getStandardTargetReflectivityData;
    end

    if nargin < 6
        CAMERA = getStandardCamera;
    end

    if nargin < 5
        alpha = 0;
    end

    if nargin < 4
        % average relative velocity in LEO is approx. 10 km/s
        v_t = 1e4;
    end

    H       = CONST.H_SOLAR;
    D       = CAMERA.D;
    Q       = CAMERA.Q;
    C       = CAMERA.C;
    chi     = CAMERA.saturation_threshold;
    P_t     = computeTransmittedPower(A_t, deltaR, alpha, D, H, TRG_DATA);
    S       = computeOpticalSignal(P_t, tau, Q, CONST);
    [L_st, W_st] = ...
            estimateStreakLengthWidth(deltaR, v_t, tau, CAMERA, CONST);
    % We check whether the streak is longer than what the detector
    % allows.
    if L_st > CAMERA.pixel_resolution * sqrt(2)
        exc = 100 * ( L_st - CAMERA.pixel_resolution * sqrt(2) ) / CAMERA.pixel_resolution * sqrt(2);
        L_st = CAMERA.pixel_resolution * sqrt(2);
        warning("WARNING: long streak, exceeds the field of view by: "+string(exc)+"%")
    end
    n_p = ceil( L_st * W_st );
    N       = computeOpticalNoise(n_p, S, tau, CAMERA);
    SNR     = S / N;
    saturation = chi * C * n_p;
    if SNR > ( saturation / N - 1 )
        %SNR = saturation / N - 1;
        %SNR = 1/ ( saturation / S - 1 );
        SNR = sqrt(saturation);
    end
end