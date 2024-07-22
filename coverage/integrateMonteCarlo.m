% integrateMonteCarlo integrates a given integrandFunction using a Monte
% Carlo method over the given sample. The result will be correct if the
% sample will be uniform over the integration volume.
%
%   integralMonteCarlo = integrateMonteCarlo(rSample, integrandFunction)
%   
%   Inputs:
%   -   rSample             : Matrix containing the integration sample. It
%                           has size nMonteCarlo x nDimensions.
%                           nMonteCarlo is the number of Monte Carlo
%                           samples, nDimensions is the number of
%                           dimensions of any given sample, i.e. the number
%                           of coordinates of any given point.
%   -   integrandFunction   : function to be integrated. Its argument is a
%                           vector of size nDimensions x 1. Its output can
%                           be a vector of any size.
%
%   Outputs:
%   -   integralMonteCarlo  : vector containing the Monte Carlo estimate
%                           for integrandFunction integral over the
%                           integration volume represented by the sample.

function integralMonteCarlo = integrateMonteCarlo(rSample, ...
                            integrandFunction)

    nMonteCarlo = size(rSample, 1);
    integralMonteCarlo = 0;
    for iMonteCarlo = 1:nMonteCarlo
        integralMonteCarlo  = integralMoneCarlo + ...
                            integrandFunction(rSample(iMonteCarlo));
    end
    integralMonteCarlo = integralMonteCarlo / nMonteCarlo;

end