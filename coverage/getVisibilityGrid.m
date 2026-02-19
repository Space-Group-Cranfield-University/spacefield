function [raMat, decMat, visibilityCountMat] = getVisibilityGrid(rObsMat, dirSun, SensorParameters, D_trg, h_trg, N_radec, R_H)
    if nargin < 7
        R_H = 6471;
    end
    if nargin < 6
        N_radec = 1e2;
    end
    if nargin < 5
        h_trg = 800;
    end
    if nargin < 4
        D_trg = 0.1;
    end
    if nargin < 3
        SensorParameters = getReducedSensorParameters;
    end
    if nargin < 2
        dirSun = [1 0 0]';
    end
    if nargin < 1
        [~, rObsMat] = initializeWalkerConstellation;
    end
    R_e = 6371;
    N = N_radec;
    raVec = linspace(0, 2*pi, 2*N);
    decVec = linspace(-pi/2, pi/2, N);
    [raMat, decMat] = meshgrid(raVec, decVec);
    visibilityCountMat = zeros(N, 2*N);
    for j = 1:(2*N)
        ra = raVec(j);
        for k = 1:N
            dec = decVec(k);
            rTrg = convertRaDecToPosition(R_e+h_trg, ra, dec);
            visibilityCountMat(k, j) = countVisibilityFOR(rTrg, rObsMat, dirSun, D_trg, SensorParameters, R_H);
        end
    end
end