function S = computeOpticalSignal(P_t, tau, Q, CONST)
    % Computes optical signal. P_t is the transmitted power, tau is the
    % exposure time.
    %
    %   S = computeOpticalSignal(P_t, tau, Q, CONST)
    %   S = computeOpticalSignal(P_t, tau, Q)
    %   S = computeOpticalSignal(P_t, tau)
    
    if nargin < 4
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 3
        Q = getStandardCamera().Q;
    end
    h = CONST.H_PLANCK;
    ni = CONST.NI_VISIBLE;
    S = P_t * tau * Q / ( h * ni );
end