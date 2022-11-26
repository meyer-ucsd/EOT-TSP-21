function [newIndexes,measurements] = getPromisingNewTargets(currentParticlesKinematic,currentParticlesExtent,currentExistences,measurements,parameters)

    % find ``free'' measurements by updating legacy potential targets only
    numMeasurements = size(measurements,2);
    numTargets = length(currentParticlesKinematic);
    measurementsCovariance = parameters.measurementVariance * eye(2);
    meanClutter = parameters.meanClutter;
    meanMeasurements = parameters.meanMeasurements;
    surveillanceRegion = parameters.surveillanceRegion;
    areaSize =  (surveillanceRegion(2,1)-surveillanceRegion(1,1)) * (surveillanceRegion(2,2)-surveillanceRegion(1,2));
    constantFactor = areaSize*(meanMeasurements/meanClutter);
    currentEpsilonExtrinsic = repmat([currentParticlesExtent]',[1,numMeasurements]);
    
    numParticles = arrayfun(@(x) length(x.particleWeights),currentParticlesKinematic);
    inputDA = zeros(numTargets,numMeasurements);
    
    for target = 1:numTargets
        lengthIndex = numParticles(target);
        targetdata=currentParticlesKinematic(:,:,target)
        targetextenddata=currentParticlesExtent(:,:,target)
        tmpMean = reshape(targetdata(1:2,:),2,lengthIndex);
        tmpCov = reshape(targetextenddata,2,2,lengthIndex) + repmat(measurementsCovariance,[1,1,lengthIndex]);
        for measurement = numMeasurements:-1:1
            likelihood1{target}(:,measurement) = constantFactor * exp(getLogWeightsFast(measurements(:,measurement),tmpMean,tmpCov));
            inputDA(target,measurement) = currentEpsilonExtrinsic(target,measurement) * currentParticlesKinematic(target).particleWeights'*likelihood1{target}(:,measurement);
        end
    end
    
    probabilitiesNew = zeros(numMeasurements,numTargets+1);
    for measurement = 1:numMeasurements
        probabilitiesNew(measurement,:) = [inputDA(:,measurement);1]/(sum(inputDA(:,measurement))+1);
    end
    probabilitiesNew = probabilitiesNew(:,end);
    
    % find central measurements in clusters of unused measurements and reorder measurement
    [newIndexes,indexesReordered] = getCentralReordered(measurements,probabilitiesNew,measurementsCovariance,parameters);
    measurements = measurements(:,indexesReordered);
    
    end
