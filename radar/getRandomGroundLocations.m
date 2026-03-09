function gndLocationMat = getRandomGroundLocations(N)
    % Requires landmask, available at: 
    % https://it.mathworks.com/matlabcentral/fileexchange/48661-landmask
    gndLocationMat = zeros(N, 2);
    k = 1;
    while k <= N
        rand_lat = rand*180 - 90;
        rand_lon = rand*360 - 180;
        if landmask(rand_lat, rand_lon, 100)
            gndLocationMat(k, :) = [rand_lon, rand_lat];
            k = k+1;
        end
    end
end