function Phi = computeMixedPhaseFunction(alpha, TRG_DATA)
    % INPUTS:   
    %   - alpha     : Phase angle 
    %   - TRG_DATA  : Target reflectivity properties
    % OUTPUTS:
    %   - Phi       : Value of the phase function at Phi(phi)
    % REQUIREMENTS:
    %   - [ 0 <= alpha <= pi ]
    %
    %   Phi = computeMixedPhaseFunction(alpha, TRG_DATA)
    %   Phi = computeMixedPhaseFunction(alpha)
    %
    % Computes the phase function of the target as observed from an optical
    % observer with phase angle phi. TRG_DATA is a data structure 
    % containing the reflectivity and mixing coefficient of the target. You
    % can initialize such data using the 
    % getStandardTargetReflectivityData() function. The model used is a 
    % mixed diffuse / specular one. The model is valid (on average) for 
    % anthropogenic space objects. If TRG_DATA is not supplied, the
    % function calls getStandardTargetReflectivityData().
    %
    % References:
    %   - 
    
    if nargin == 1
        TRG_DATA = getStandardTargetReflectivityData;
    end
    % assert((alpha >= 0 && alpha <= pi), "Phase angle should be between 0 deg and 180 deg")
    rho     = TRG_DATA.rho_B;
    k       = TRG_DATA.k_m;
    K_1     = k * 2 / ( 3 * pi^2 );
    K_2     = ( 1 - k ) / 4;
    Phi     = rho * ...
            ( K_1 * ( ( pi - alpha ) .* cos(alpha) + sin(alpha) ) + K_2 );
end