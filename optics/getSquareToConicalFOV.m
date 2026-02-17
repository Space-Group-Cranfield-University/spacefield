function FOV_eq = getSquareToConicalFOV(FOV)
    % Computes the equivalent conical FOV from the original FOV for a
    % square detector.
    FOV_eq = sqrt(4/pi) * FOV;
end