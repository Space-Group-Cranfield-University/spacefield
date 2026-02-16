function N = estimateStreaksWithinFOV(tau, D_1, D_2, FOV, v_t, deltaRmax)
    % Estimates the number of streaks falling within the FOV over the
    % observation exposure time tau for objects of sizes between D_1 and 
    % D_2. Valid for sizes > 1 mm and < 10 m. Assumes average
    % target-observer relative velocity parallel to the FOV plane is 9.5
    % km/s if not supplied.
    if nargin < 6
        deltaRmax = 42000000;
    end
    if nargin < 5
        v_t = 9500;
    end
    if nargin < 4
        FOV = getStandardCamera().FOV(1);
    end
    if nargin < 2
        D_1 = 1e-2;
        D_2 = 10;
    end
    
    V_lim = estimateLimitingMagnitude(tau);
    K = ( 26.74 + V_lim + 2.5 * log10(computeMixedPhaseFunction(0)) ) / 5;
    D_sharp = 10^(-K + log10(deltaRmax));
    % disp('V_lim [m]: '+string(V_lim)+', D_sharp [m]: '+string(D_sharp))
    phi_1 = 10^( - 6 - log10(D_1) );
    phi_2 = 10^( - 6 - log10(D_2) );
    phi_sharp = 10^( - 6 - log10(D_sharp) );
    if D_sharp >= D_2
        N   = pi * (FOV / 2) * ...
            (   10^( -6 + 2 * K ) * ( D_2 - D_1 ) ...
                - ( v_t^2 .* tau.^2 ./ FOV^2 ) * ( phi_1 - phi_2 ) ...
            );
    elseif D_sharp > D_1
        N   = pi * (FOV / 2) * ...
            (   10^( -6 + 2 * K ) * ( D_sharp - D_1 ) ...
                - ( v_t^2 .* tau.^2 ./ FOV^2 ) * ( phi_1 - phi_2 ) ...
                + deltaRmax^2 * (phi_sharp - phi_2) ...
            );
    else
        N   = pi * (FOV / 2) * ...
            (   - ( v_t^2 .* tau.^2 ./ FOV^2 ) * ( phi_1 - phi_2 ) ...
                + deltaRmax^2 * (phi_sharp - phi_2) ...
            );
    end
    N = N/(86400 * 365) .* tau;
end