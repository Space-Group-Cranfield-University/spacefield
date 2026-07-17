function DCM = getRotationMatrixLVLHToECI(x)
    r = x(1:3);
    v = x(4:6);
    h = cross(r, v);
    e_x = r / norm(r);
    e_z = h / norm(h);
    e_y = cross(e_z, e_x);
    DCM = [e_x, e_y, e_z];
end