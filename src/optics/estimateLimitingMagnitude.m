function V_lim = estimateLimitingMagnitude(CAMERA, tau, SNR_threshold, N, N_s)
    % This function estimates the limiting magnitude of a sensor through
    % Monte Carlo sampling of target parameters. Considers the SNR
    % threhsold equal to 5.
    %
    % N     : Number of Monte Carlo samples (1e4 is usually fast and 
    %       accurate)
    % N_s   : Number of subsamples used to estimate limiting magnitude.
    %       Needs to be <= N!
    if nargin < 1
        CAMERA = getStandardCamera();
    end
    if nargin < 5
        N = 1e4;
        N_s = 20;
    end
    if nargin < 3
        SNR_threshold = CAMERA.SNR_threshold;
    end
    if nargin < 2
        tau = CAMERA.tau;
    end
    deltaRvecExp = [5 7];
    D_tVecExp = [-2 -1];
    v_tVecExp = [2 3];
    sampleExp = @(p) 10.^(p(1) + rand*(p(2)-p(1)));
    
    for k = 1:N
        deltaR = sampleExp(deltaRvecExp);
        v_t = sampleExp(v_tVecExp);
        D_t = sampleExp(D_tVecExp);
        A_t = D_t^2;
        alpha = rand*pi;
        SNR(k)  = computeOpticalSNR(tau, deltaR, A_t, v_t, alpha, CAMERA);
        V(k)    = computeApparentMagnitude(deltaR, A_t, alpha);
        dist(k) = abs(SNR(k) - SNR_threshold);
    end
    [~, index] = sort(dist);
    V = V(index);
    V_lim = 0;
    for k = 1:N_s
        V_lim = V_lim + V(k);
    end
    V_lim = V_lim / N_s;
end