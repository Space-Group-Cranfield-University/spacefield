function [p, rsq] = linearRegression(xData, yData)
    p = polyfit(xData, yData, 1);
    yfit = polyval(p, xData);
    yresid = yData - yfit;
    SSresid = sum(yresid.^2);
    SStotal = sum(length(yData)-1) * var(yData);
    rsq = 1 - SSresid/SStotal;
end