function rTrgAugmentedMat = sampleTargetPopulationFromDistribution(nTrg, fileNames, dataFolder, R_E)
    % Samples nTrg targets from a given altitude / declination / size
    % distribution. The function is not very modular because randsample()
    % works better in vectorized form.
    %
    % Requires: Statistics and Machine Learrning Toolbox

    if nargin < 3
        dataFolder = 'data';
    end
    if nargin < 4
        R_E = 6371;
    end
    
    AltitudeData = importAltitudeData(fileNames(1), dataFolder);
    DeclinationData = importDeclinationData(fileNames(2), dataFolder);
    SizeData = importSizeData(fileNames(3), dataFolder);
    AltitudeData{:,2} = AltitudeData{:,2}/sum(AltitudeData{:,2});
    DeclinationData{:,2} = DeclinationData{:,2}/sum(DeclinationData{:,2});
    SizeData{:,2} = SizeData{:,2}/(sum(SizeData{:,2})); 
    sizeSample = randsample(SizeData{:,1}, nTrg, true, SizeData{:,2});
    altitudeSample = randsample(AltitudeData{:,1}, nTrg, true, AltitudeData{:,2});
    declinationSample = deg2rad(randsample(DeclinationData{:,1}, nTrg, true, DeclinationData{:,2}));
    sample = zeros(nTrg, 4);
    for k = 1:nTrg
        raSampleTemp = rand()*2*pi;
        sizeSampleTemp = sizeSample(k);
        altitudeSampleTemp = altitudeSample(k);
        declinationSampleTemp = declinationSample(k);
        radiusTemp = altitudeSampleTemp + R_E;
        x = radiusTemp * cos(raSampleTemp) * cos(declinationSampleTemp);
        y = radiusTemp * sin(raSampleTemp) * cos(declinationSampleTemp);
        z = radiusTemp * sin(declinationSampleTemp);
        sampleTemp = [x, y, z, sizeSampleTemp];
        sample(k, :) = sampleTemp;
    end
    rTrgAugmentedMat = sample; % Store the sampled targets in the output variable
end