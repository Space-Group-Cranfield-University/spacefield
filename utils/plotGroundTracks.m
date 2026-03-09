function [] = plotGroundTracks(raDecMat, color, lineWidth, Z)
    if nargin < 4
        Z = 100;
    end
    if nargin < 3
        lineWidth = 1;
    end
    if nargin < 2
        color = "red";
    end
    hold on
    nRaDec = size(raDecMat,1);
    raDecPlot = raDecMat(1,2:3);
    for k = 2:nRaDec
        if abs(raDecMat(k,2) - raDecMat(k-1,2)) > pi
            plot3(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), Z*ones(size(raDecPlot, 1),1), "Color", color, "LineWidth", lineWidth)
            raDecPlot = raDecMat(k,2:3);
        else
            raDecPlot = [raDecPlot; raDecMat(k,2:3)];
        end
    end
    plot3(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), Z*ones(size(raDecPlot, 1),1), "Color", color, "LineWidth", lineWidth)
end