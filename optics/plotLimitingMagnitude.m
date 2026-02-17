function [] = plotLimitingMagnitude(tau, SNR_threshold, CAMERA, N)
    if nargin < 4
        N = 1e3;
    end
    if nargin < 3
        CAMERA = getStandardCamera();
    end
    if nargin < 2
        SNR_threshold = 5;
    end
    if nargin < 1
        tau = 0.1;
    end
    sampleExp = @(p) 10.^(p(1) + rand*(p(2)-p(1)));
    deltaRvecExp = [3 7];
    D_tVecExp = [-2 1];
    v_tVecExp = [2 3];
    for k = 1:N
        deltaR = sampleExp(deltaRvecExp);
        v_t = sampleExp(v_tVecExp);
        D_t = sampleExp(D_tVecExp);
        A_t = D_t^2;
        alpha = rand*pi;
        SNR(k)  = computeOpticalSNR(tau, deltaR, A_t, v_t, alpha, CAMERA);
        V(k)    = computeApparentMagnitude(deltaR, A_t, alpha);
    end
    minSNR = 10*log10(SNR_threshold);
    f = figure;
    hold on
    scatter(V, 10*log10(SNR), "filled")
    plot([-15 22], [minSNR minSNR], "red")
    text(-10, 10-10*log10(5)+10*log10(SNR_threshold),"SNR = "+string(SNR_threshold),'Color','red', 'LineWidth', 2, 'FontSize', 14);
    xlabel("apparent magnitude", "FontSize", 14)
    ylabel("SNR [dB]", "FontSize", 14)
    %legend(legendVec, "FontSize", 14, 'Interpreter','latex')
    set(gca,'fontsize',14)
end