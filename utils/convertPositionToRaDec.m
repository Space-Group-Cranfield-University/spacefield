function raDec = convertPositionToRaDec(r, twoPiFlag)
    if nargin < 2
        twoPiFlag = 0;
    end
    x = r(1);
    y = r(2);
    z = r(3);
    rNorm = norm(r);
    dec = asin( z / rNorm );
    ra = atan2( y, x );
    if twoPiFlag && ra < 0
        ra = 2*pi - abs(ra);
    end
    raDec = [rNorm; ra ; dec];
end