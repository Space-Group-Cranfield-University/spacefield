function raDec = convertPositionToRaDec(r)
    x = r(1);
    y = r(2);
    z = r(3);
    rNorm = norm(r);
    dec = asin( z / rNorm );
    ra = atan2( y, x );
    if ra < 0
        ra = 2*pi - abs(ra);
    end
    raDec = [ra ; dec; rNorm];
end