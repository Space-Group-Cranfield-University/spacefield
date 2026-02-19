function rTrg = convertRaDecToPosition(R, ra, dec)
    rTrg = R * ones(3, 1);
    rTrg(1) = rTrg(1) * cos(dec) * cos(ra);
    rTrg(2) = rTrg(2) * cos(dec) * sin(ra);
    rTrg(3) = rTrg(3) * sin(dec);
end