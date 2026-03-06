function rTrg = convertRaDecToPosition(raDec, ra, dec)
    if nargin == 1
        R = raDec(1);
        ra = raDec(2);
        dec = raDec(3);
    else
        R = raDec;
    end
    rTrg = R * ones(3, 1);
    rTrg(1) = rTrg(1) * cos(dec) * cos(ra);
    rTrg(2) = rTrg(2) * cos(dec) * sin(ra);
    rTrg(3) = rTrg(3) * sin(dec);
end