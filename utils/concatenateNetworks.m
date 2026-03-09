function NETWORK = concatenateNetworks(GND, OPT)
    OPT = rmfield(OPT, "a");
    OPT = rmfield(OPT, "e");
    OPT = rmfield(OPT, "in");
    OPT = rmfield(OPT, "raan");
    OPT = rmfield(OPT, "aop");
    OPT = rmfield(OPT, "v");
    OPT = rmfield(OPT, "kep");
    OPT = rmfield(OPT, "T");
    GND = rmfield(GND, "t_0");
    GND = rmfield(GND, "lonLat");
    GND = rmfield(GND, "r0_sph");
    GND = rmfield(GND, "r0_ECEF");
    GND = rmfield(GND, "r0_ECI");
    GND = rmfield(GND, "v0");
    NETWORK = [GND, OPT];
end