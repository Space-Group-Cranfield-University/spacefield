function DCM = getRotationMatrix_X(theta)
    c = cos(theta);
    s = sin(theta);
    DCM = [1 0 0;
        0 c -s;
        0 s c];
end