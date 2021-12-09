% Florian Meyer, 2020

function [weights,updatedExistence] = getWeightsUnknown(logWeights,oldExistence,skipIndex)

if(skipIndex)
    logWeights(:,skipIndex) = zeros(size(logWeights,1),1);
end

logWeights = sum(logWeights,2);

aliveUpdate = mean(exp(logWeights),1);
if(isinf(aliveUpdate))
    updatedExistence = 1;
else
    alive = oldExistence * aliveUpdate;
    dead = (1 - oldExistence);
    updatedExistence = alive / (dead + alive);
end

weights = exp(logWeights - max(logWeights));
weights = 1/sum(weights,1) * weights;

end