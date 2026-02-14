function P_t = computeTransmittedPower(A_t, deltaR, alpha, D, H, TRG_DATA)
    %   Computes power transmitted from a space target to an optical
    %   sensor. Uses a mixed diffuse / specular model for reflection. A_t
    %   is the target cross-section area, deltaR is the relative distance
    %   in m, alpha is the phase angle.
    %
    %     P_t = computeTransmittedPower(A_t, deltaR, alpha, D, H, TRG_DATA)
    %     P_t = computeTransmittedPower(A_t, deltaR, alpha, D, H)
    %     P_t = computeTransmittedPower(A_t, deltaR, alpha, D)
    %     P_t = computeTransmittedPower(A_t, deltaR, alpha)

    if nargin < 6
        Phi = computeMixedPhaseFunction(alpha);
    else
        Phi = computeMixedPhaseFunction(alpha, TRG_DATA);
    end

    if nargin < 5
        H = initialiseAstronomicalConstants().H_SOLAR;
    end

    if nargin < 4
        D = getStandardCamera().D;
    end

    P_t = H * A_t * Phi * pi * D^2 / ( 4 * deltaR^2 );

end