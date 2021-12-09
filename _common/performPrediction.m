% Florian Meyer, 2020

function [newParticles, newExistences, newExtents] = performPrediction( oldParticles, oldExistences, oldExtents, scanTime, parameters )
[~,numParticles,numTargets] = size(oldParticles);
drivingNoiseVariance = parameters.accelerationDeviation^2;
survivalProbability = parameters.survivalProbability;
degreeFreedomPrediction = parameters.degreeFreedomPrediction;

[A, W] = getTransitionMatrices(scanTime);
newParticles = oldParticles;
newExistences = oldExistences;
newExtents = oldExtents;

for target = 1:numTargets
    oldExtents(:,:,:,target) = oldExtents(:,:,:,target)/degreeFreedomPrediction;
    newExtents(:,:,:,target) = wishrndFastVector(oldExtents(:,:,:,target),degreeFreedomPrediction,numParticles);
end

for target = 1:numTargets
    newParticles(:,:,target) = A*oldParticles(:,:,target) + W*sqrt(drivingNoiseVariance)*randn(2,numParticles);
    newExistences(target) = survivalProbability*oldExistences(target);
end

end

function [A, W] = getTransitionMatrices( scanTime )
A = diag(ones(4,1));
A(1,3) = scanTime;
A(2,4) = scanTime;

W = zeros(4,2);
W(1,1) = 0.5*scanTime^2;
W(2,2) = 0.5*scanTime^2;
W(3,1) = scanTime;
W(4,2) = scanTime;
end