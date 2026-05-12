function DCM = getRotationMatrix_Y(theta)
    c = cos(theta);
    s = sin(theta);
    DCM = [c 0 s;
        0 1 0;
        -s 0 c];
end