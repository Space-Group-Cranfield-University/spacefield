function [L_st, W_st] = estimateStreakLengthWidth(deltaR, v_t, tau, CAMERA, CONST)
    % Estimates the number of pixels covered by the streak produced by a
    % space target over a camera field-of-view. deltaR is the relative
    % distance in m, v_t is the apparent tangential velocity in m/s, and 
    % tau is the exposure time in s. The outpur lengths are in pixel units.
    %
    % [L_st, W_st] = estimateStreakLengthWidth(deltaR, v_t, tau, CAMERA, CONST)
    % [L_st, W_st] = estimateStreakLengthWidth(deltaR, v_t, tau, CAMERA)
    % [L_st, W_st] = estimateStreakLengthWidth(deltaR, v_t, tau)
    
    if nargin < 5
        lambda = initialiseAstronomicalConstants().LAMBDA_VISIBLE;
    else
        lambda = CONST.LAMBDA_VISIBLE;
    end

    if nargin < 4
        CAMERA = getStandardCamera;
    end
    f = CAMERA.f;
    delta_p = CAMERA.pixel_size;
    F = CAMERA.F;
    L_st = ceil( f * v_t * tau / ( delta_p * deltaR ) );
    W_st = ceil( 2.44 * lambda * F / delta_p );
end