function D = getStandardScaleMatrix()
    D = blkdiag(eye(3)*(1/7500), eye(3)*(1/7.5));
end