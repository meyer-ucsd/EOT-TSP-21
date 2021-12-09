% Florian Meyer, 2020

function [ estimatedTracks, estimatedExtents ] = eotEllipticalShape( measurementsCell, parameters )
numParticles = parameters.numParticles;
meanClutter = parameters.meanClutter;
meanMeasurements = parameters.meanMeasurements;
scanTime = parameters.scanTime;
detectionThreshold = parameters.detectionThreshold;
thresholdPruning = parameters.thresholdPruning;
numOuterIterations = parameters.numOuterIterations;
meanBirths = parameters.meanBirths;
priorVelocityCovariance = parameters.priorVelocityCovariance;
surveillanceRegion = parameters.surveillanceRegion;

areaSize =  (surveillanceRegion(2,1)-surveillanceRegion(1,1)) * (surveillanceRegion(2,2)-surveillanceRegion(1,2));
measurementsCovariance = parameters.measurementVariance * eye(2);

priorExtent1 = parameters.priorExtent1;
priorExtent2 = parameters.priorExtent2;
meanExtentPrior = (parameters.priorExtent1/(parameters.priorExtent2-3))^2;
totalCovariance = meanExtentPrior^2+measurementsCovariance;

[numSteps, ~] = size(measurementsCell);
constantFactor = areaSize*(meanMeasurements/meanClutter);
uniformWeight = log(1/areaSize);


estimates = cell(numSteps,1);
currentLabels = zeros(2,0);
currentParticlesKinematic = zeros(4,numParticles,0);
currentExistences = zeros(0,1);
currentParticlesExtent = zeros(2,2,numParticles,0);

for step = 1:numSteps
    step
    
    % load current measurements
    measurements = measurementsCell{step};
    numMeasurements = size(measurements,2);
    
    
    % perform prediction step
    [currentParticlesKinematic,currentExistences,currentParticlesExtent] = performPrediction(currentParticlesKinematic,currentExistences,currentParticlesExtent,scanTime,parameters);    
    currentAlive = currentExistences*exp(-meanMeasurements);
    currentDead = (1-currentExistences);
    currentExistences = currentAlive./(currentDead+currentAlive);
    numTargets = size(currentParticlesKinematic,3);
    numLegacy = numTargets;
    
    
    % get indexes of promising new objects 
    [newIndexes,measurements] = getPromisingNewTargets(currentParticlesKinematic,currentParticlesExtent,currentExistences,measurements,parameters);
    numNew = size(newIndexes,1);
    currentLabels = cat(2,currentLabels,[step*ones(1,numNew);newIndexes']);
    
    
    % initialize belief propagation (BP) message passing
    newExistences = repmat(meanBirths * exp(-meanMeasurements)/(meanBirths * exp(-meanMeasurements) + 1),[numNew,1]);
    newParticlesKinematic = zeros(4,numParticles,numNew);
    newParticlesExtent = zeros(2,2,numParticles,numNew);
    newWeights = zeros(numParticles,numNew);
    for target = 1:numNew
        proposalMean = measurements(:,newIndexes(target));
        proposalCovariance = 2 * totalCovariance; % strech covariance matrix to make proposal distribution heavier-tailed then target distribution
        
        newParticlesKinematic(1:2,:,target) = proposalMean + sqrtm(proposalCovariance) * randn(2,numParticles);
        newWeights(:,target) = uniformWeight - log(mvnpdf(newParticlesKinematic(1:2,:,target)', proposalMean', proposalCovariance));
        
        newParticlesExtent(:,:,:,target) = iwishrndFastVector(priorExtent1,priorExtent2,numParticles);
    end
    
    currentExistences = cat(1,currentExistences,newExistences);
    currentExistencesExtrinsic = repmat(currentExistences,[1,numMeasurements]);
    
    currentParticlesKinematic = cat(3,currentParticlesKinematic,newParticlesKinematic);
    currentParticlesExtent = cat(4,currentParticlesExtent,newParticlesExtent);
    
    weightsExtrinsic = nan(numParticles,numMeasurements,numLegacy);
    weightsExtrinsicNew = nan(numParticles,numMeasurements,size(newIndexes,1));
    
    likelihood1 = zeros(numParticles,numMeasurements,numTargets);
    likelihoodNew1 = nan(numParticles,numMeasurements,size(newIndexes,1));
    for outer = 1:numOuterIterations
        
        % perform one BP message passing iteration for each measurement
        outputDA = cell(numMeasurements,1);
        targetIndexes = cell(numMeasurements,1);
        for measurement = numMeasurements:-1:1
            inputDA = ones(2,numLegacy);
            
            for target = 1:numLegacy
                
                if(outer == 1)
                    likelihood1(:,measurement,target) = constantFactor * exp(getLogWeightsFast(measurements(:,measurement),currentParticlesKinematic(1:2,:,target),getSquare2Fast(currentParticlesExtent(:,:,:,target)) + repmat(measurementsCovariance,[1,1,numParticles])));
                    inputDA(2,target) = currentExistencesExtrinsic(target,measurement) * mean(likelihood1(:,measurement,target),1);
                else
                    inputDA(2,target) = currentExistencesExtrinsic(target,measurement) * (weightsExtrinsic(:,measurement,target)'*likelihood1(:,measurement,target));
                end
                
                inputDA(1,target) = 1;
            end
            
            targetIndex = numLegacy;
            targetIndexesCurrent = nan(numLegacy,1);
            
            % only new targets with index >= measurement index are connected to measurement
            for target = numMeasurements:-1:measurement
                
                if(any(target==newIndexes))
                    targetIndex = targetIndex + 1;
                    targetIndexesCurrent = [targetIndexesCurrent;target];
                    
                    if(outer == 1)
                        weights = exp(newWeights(:,targetIndex-numLegacy));
                        weights = (weights/sum(weights,1))';
                        likelihoodNew1(:,measurement,targetIndex-numLegacy) = constantFactor * exp(getLogWeightsFast(measurements(:,measurement),currentParticlesKinematic(1:2,:,targetIndex),getSquare2Fast(currentParticlesExtent(:,:,:,targetIndex)) + repmat(measurementsCovariance,[1,1,numParticles])));
                        inputDA(2,targetIndex) = currentExistencesExtrinsic(targetIndex,measurement) * (weights*likelihoodNew1(:,measurement,targetIndex-numLegacy));
                    else
                        inputDA(2,targetIndex) = currentExistencesExtrinsic(targetIndex,measurement) *  (weightsExtrinsicNew(:,measurement,targetIndex-numLegacy)'*likelihoodNew1(:,measurement,targetIndex-numLegacy));
                    end
                    inputDA(1,targetIndex) = 1;
                    
                    if(target == measurement)
                        inputDA(1,targetIndex) = 1 - currentExistencesExtrinsic(targetIndex,measurement);
                    end
                end
            end
           
         targetIndexes{measurement} = targetIndexesCurrent;   
         outputDA{measurement} = dataAssociationBP(inputDA);   
         
        end
        
        
        % perform update step for legacy targets
        for target = 1:numLegacy
            weights = zeros(size(currentParticlesKinematic,2),numMeasurements);
            for measurement = 1:numMeasurements
                currentWeights = 1 + likelihood1(:,measurement,target) * outputDA{measurement}(1,target);
                currentWeights = log(currentWeights);
                weights(:,measurement) = currentWeights;
            end
            
            % calculate extrinsic information for legacy targets (at all except last iteration) and belief (at last iteration)
            if(outer ~= numOuterIterations)
                for measurement = 1:numMeasurements
                    [weightsExtrinsic(:,measurement,target),currentExistencesExtrinsic(target,measurement)] = getWeightsUnknown(weights,currentExistences(target),measurement);
                end
            else
                [currentParticlesKinematic(:,:,target),currentParticlesExtent(:,:,:,target),currentExistences(target)] = updateParticles(currentParticlesKinematic(:,:,target),currentParticlesExtent(:,:,:,target),currentExistences(target),weights,parameters);
            end
        end
        
        % perform update step for new targets
        targetIndex = numLegacy;
        for target = numMeasurements:-1:1
            if(any(target == newIndexes))
                
                targetIndex = targetIndex + 1;
                weights = zeros(size(currentParticlesKinematic,2),numMeasurements+1);
                weights(:,numMeasurements+1) = newWeights(:,targetIndex-numLegacy);
                for measurement = 1:target
                    
                    outputTmpDA = outputDA{measurement}(1,targetIndexes{measurement}==target);
                    
                    if(~isinf(outputTmpDA))
                        currentWeights = likelihoodNew1(:,measurement,targetIndex-numLegacy) * outputTmpDA;
                    else
                        currentWeights = likelihoodNew1(:,measurement,targetIndex-numLegacy);
                    end
                    
                    if(measurement ~= target)
                        currentWeights = currentWeights + 1;
                    end
                    currentWeights = log(currentWeights);
                    weights(:,measurement) = currentWeights;
                end
                
                % calculate extrinsic information for new targets (at all except last iteration) or belief (at last iteration)
                if(outer ~= numOuterIterations)
                    for measurement = 1:target
                        [weightsExtrinsicNew(:,measurement,targetIndex-numLegacy),currentExistencesExtrinsic(targetIndex,measurement)] = getWeightsUnknown(weights,currentExistences(targetIndex),measurement);
                    end
                else
                    [currentParticlesKinematic(1:2,:,targetIndex),currentParticlesExtent(:,:,:,targetIndex),currentExistences(targetIndex)] = updateParticles(currentParticlesKinematic(1:2,:,targetIndex),currentParticlesExtent(:,:,:,targetIndex),currentExistences(targetIndex),weights,parameters);
                    currentParticlesKinematic(3:4,:,targetIndex) = mvnrnd([0;0],priorVelocityCovariance,numParticles)';
                end
            end
        end
    end
    
    % perform pruning
    numTargets = size(currentParticlesKinematic,3);
    isRedundant = false(numTargets,1);
    for target = 1:numTargets
        if(currentExistences(target) < thresholdPruning)
            isRedundant(target) = true;
        end
    end
    currentParticlesKinematic = currentParticlesKinematic(:,:,~isRedundant);
    currentParticlesExtent = currentParticlesExtent(:,:,:,~isRedundant);
    currentLabels = currentLabels(:,~isRedundant);
    currentExistences = currentExistences(~isRedundant);
    
    
    % perform estimation
    numTargets = size(currentParticlesKinematic,3);
    detectedTargets = 0;
    for target = 1:numTargets
        if(currentExistences(target) > detectionThreshold)
            detectedTargets = detectedTargets + 1;
            estimates{step}.state(:,detectedTargets) = mean(currentParticlesKinematic(:,:,target),2);
            estimates{step}.extent(:,:,detectedTargets) = mean(currentParticlesExtent(:,:,:,target),3);
            estimates{step}.label(:,detectedTargets) = currentLabels(:,target);
        end
    end
end

[estimatedTracks,estimatedExtents] = trackFormation(estimates, parameters);
end