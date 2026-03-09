function r_2 = applyRodriguezFormula(r_1, h, theta)
    % https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    r_2 = r_1 * cos(theta) + cross(h, r_1) * sin(theta) + h * (h' * r_1) * (1 - cos(theta));
end