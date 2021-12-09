function [targetTracks, targetExtents] = generateTracksUnknown(parameters,startStates,extentMatrixes,appearanceFromTo,numSteps)
accelerationDeviation = parameters.accelerationDeviation;
scanTime = parameters.scanTime;
numTargets = size(startStates,2);

[A, W] = getTransitionMatrices(scanTime);
targetTracks = nan(4,numSteps,numTargets);
targetExtents = nan(2,2,numSteps,numTargets);
for target = 1:numTargets
    tmp = startStates(:,target);
    for step = 1:numSteps
        tmp = A*tmp + W*accelerationDeviation*randn(2,1);
        if( (step >= appearanceFromTo(1,target)) && (step <= appearanceFromTo(2,target)) )
            targetTracks(:,step,target) = tmp;
            targetExtents(:,:,step,target) = extentMatrixes(:,:,target);
        end
    end
end

end

