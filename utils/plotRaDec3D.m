function [] = plotRaDec3D(raDecMat, color, lineWidth, Z)
    hold on
    nRaDec = size(raDecMat,1);
    raDecPlot = raDecMat(1,:);
    for k = 2:nRaDec
        if abs(raDecMat(k,1) - raDecMat(k-1,1)) > pi
            plot3(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), Z*ones(size(raDecPlot, 1),1), "Color", color, "LineWidth", lineWidth)
            raDecPlot = raDecMat(k,:);
        else
            raDecPlot = [raDecPlot; raDecMat(k,:)];
        end
    end
    plot3(rad2deg(raDecPlot(:,1)), rad2deg(raDecPlot(:,2)), Z*ones(size(raDecPlot, 1),1), "Color", color, "LineWidth", lineWidth)
end