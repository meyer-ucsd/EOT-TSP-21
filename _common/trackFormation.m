% Florian Meyer, 2020

function [tracks,extents] = trackFormation(estimates, parameters)
numSteps = length(estimates);
minimumTrackLength = parameters.minimumTrackLength;

% get all labels
labels = zeros(2,0);
for step = 1:numSteps
    if(~isempty(estimates{step}))
        currentLabels = [estimates{step}.label];
        labels = [labels,currentLabels];
    end
end
labels = unique(labels','rows')';

% write estimated states in the correct track
[~,numTracks] = size(labels);
tracks = nan(4,numSteps,numTracks);
extents = nan(2,2,numSteps,numTracks);
for step = 1:numSteps
    if(~isempty(estimates{step}))
        currentLabels = [estimates{step}.label];
        currentStates = [estimates{step}.state];
        currentExtents = [estimates{step}.extent];
        [~, numEstimates] = size(currentLabels);
        for estimate = 1:numEstimates
            indexes =  compareVectorWithMatrix( labels, currentLabels(:,estimate) );
            tracks(:,step,indexes) = currentStates(:,estimate);
            extents(:,:,step,indexes) = currentExtents(:,:,estimate);
        end
    end
end

% remove too short tracks
indexes = false(numTracks,1);
for track = 1:numTracks
    if(sum(~isnan(tracks(1,:,track))) >= minimumTrackLength)
        indexes(track) = true;
    end
end
tracks = tracks(:,:,indexes);
extents = extents(:,:,:,indexes);

end

function [ indexes ] = compareVectorWithMatrix( matrix, vector )
indexes = [];

[numRows, numColums] = size(matrix);
[numRowsTmp, numColumsTmp] = size(vector);

if(numRows==numRowsTmp && numColumsTmp == 1)
    checkmatrix = (matrix==repmat(vector,[1,numColums]));
    indexes = (sum(checkmatrix)==numRows);
end

end