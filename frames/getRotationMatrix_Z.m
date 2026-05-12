function DCM = getRotationMatrix_Z(theta)
    c = cos(theta);
    s = sin(theta);
    DCM = [c -s 0;
        s c 0;
        0 0 1];
end