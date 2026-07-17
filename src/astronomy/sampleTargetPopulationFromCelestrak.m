function TRG = sampleTargetPopulationFromCelestrak(N_t, DATA, sizeFileName, dataFolder, MU, R)
    if nargin < 5
        MU = initialiseAstronomicalConstants().MU_E;
        R = initializeAstronomicalConstants().R_E;
    end
    if nargin < 4
        dataFolder = "data";
    end
    if isnumeric(sizeFileName)
        sizeSample = sizeFileName * ones(1, N_t);
    else
        SizeData = importSizeData(sizeFileName, dataFolder);
        SizeData{:,2} = SizeData{:,2}/(sum(SizeData{:,2}));
        sizeSample = randsample(SizeData{:,1}, N_t, true, SizeData{:,2});
    end
    TRG = struct('a', cell(1, N_t), 'e', cell(1, N_t),'in', cell(1, N_t),...
        'raan',cell(1, N_t),'aop', cell(1, N_t),'v', cell(1, N_t),...
        'kep',cell(1, N_t),'x0', cell(1, N_t));
    for k = 1:N_t
        index = randi(size(DATA, 1));
        a   = R + ( DATA.APOGEE(index) + DATA.PERIGEE(index) ) / 2;
        e   = (DATA.APOGEE(index) - DATA.PERIGEE(index)) / ...
            (2 * R + DATA.APOGEE(index) + DATA.PERIGEE(index));
        in  = deg2rad(DATA.INCLINATION(index));
        TRG(k).a = a;
        TRG(k).e = e;
        TRG(k).in = in;
        TRG(k).raan = rand * 2 * pi;
        TRG(k).aop = rand * 2 * pi;
        TRG(k).v = rand * 2 * pi;
        TRG(k).kep  = [TRG(k).a, TRG(k).e, TRG(k).in, TRG(k).raan, ...
                    TRG(k).aop, TRG(k).v]';
        TRG(k).x0 = convertKepToCart(TRG(k).kep, MU);
        TRG(k).D_t = sizeSample(k);
    end
end