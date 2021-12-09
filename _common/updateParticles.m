% Florian Meyer, 2020

function [updatedParticlesKinematic,updatedParticlesExtent,updatedExistence] = updateParticles(oldParticlesKinematic,oldParticlesExtent,oldExistence,logWeights,parameters)
numParticles = parameters.numParticles;
regularizationDeviation = parameters.regularizationDeviation;


logWeights = sum(logWeights,2);

aliveUpdate = mean(exp(logWeights),1);
if(isinf(aliveUpdate))
    updatedExistence = 1;
else
    alive = oldExistence*aliveUpdate;
    dead = (1-oldExistence);
    updatedExistence = alive/(dead+alive);
end

if(updatedExistence ~= 0)
    logWeights = logWeights-max(logWeights);
    weights = exp(logWeights);
    weightsNormalized = 1/sum(weights)*weights;
    
    indexes = resampleSystematic(weightsNormalized,numParticles);
    updatedParticlesKinematic = oldParticlesKinematic(:,indexes);
    updatedParticlesExtent = oldParticlesExtent(:,:,indexes);
    
    updatedParticlesKinematic(1:2,:) = updatedParticlesKinematic(1:2,:) + regularizationDeviation * randn(2,numParticles);
else
    updatedParticlesKinematic = nan(size(oldParticlesKinematic));
    updatedParticlesExtent = nan(size(oldParticlesExtent));
end

end