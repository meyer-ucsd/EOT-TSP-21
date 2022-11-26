function [ newIndexes, measurements ] = getPromisingNewTargets( currentParticlesKinematicTmp, currentParticlesExtentTmp, currentExistencesTmp, measurements, parameters )
numMeasurements = size( measurements, 2 );
numParticles = size( currentParticlesKinematicTmp, 2 );
measurementsCovariance = parameters.measurementVariance * eye( 2 );
surveillanceRegion = parameters.surveillanceRegion;
areaSize = ( surveillanceRegion( 2, 1 ) - surveillanceRegion( 1, 1 ) ) * ( surveillanceRegion( 2, 2 ) - surveillanceRegion( 1, 2 ) );
meanMeasurements = parameters.meanMeasurements;
meanClutter = parameters.meanClutter;
constantFactor = areaSize * ( meanMeasurements / meanClutter );

probabilitiesNew = ones( numMeasurements, 1 );
for measurement = 1:numMeasurements
    numTargets = size( currentParticlesKinematicTmp, 3 );

    inputDA = ones( 2, numTargets );
    likelihoods = zeros( numParticles, numTargets );
    for target = 1:numTargets
        likelihoods( :, target ) = constantFactor .* exp( getLogWeightsFast( measurements( :, measurement ), currentParticlesKinematicTmp( 1:2, :, target ), getSquare2Fast( currentParticlesExtentTmp( :, :, :, target ) ) + repmat( measurementsCovariance, [ 1, 1, numParticles ] ) ) );
        inputDA( 2, target ) = mean( likelihoods( :, target ), 1 );
        inputDA( :, target ) = currentExistencesTmp( target ) * inputDA( :, target ) + ( 1 - currentExistencesTmp( target ) ) * [ 1;0 ];
    end

    inputDA = inputDA( 2, : ) ./ inputDA( 1, : );
    sumInputDA = 1 + sum( inputDA, 2 );
    outputDA = 1 ./ ( repmat( sumInputDA, [ 1, numTargets ] ) - inputDA );
    probabilitiesNew( measurement ) = 1 ./ sumInputDA;

    if ( measurement == numMeasurements )
        break ;
    end

    for target = 1:numTargets
        logWeights = log( ones( numParticles, 1 ) + likelihoods( :, target ) * outputDA( 1, target ) );
        [ currentParticlesKinematicTmp( :, :, target ), currentParticlesExtentTmp( :, :, :, target ), currentExistencesTmp( target ) ] = updateParticles( currentParticlesKinematicTmp( :, :, target ), currentParticlesExtentTmp( :, :, :, target ), currentExistencesTmp( target ), logWeights, parameters );
    end
end

[ newIndexes, indexesReordered ] = getCentralReordered( measurements, probabilitiesNew, parameters );

measurements = measurements( :, indexesReordered );

end

function [ centralIndexes, indexesReordered ] = getCentralReordered( measurements, probabilitiesNew, parameters )
parameters.freeThreshold = .85;
parameters.clusterThreshold = 0.9;
parameters.minClusterElements = 2;

threshold = parameters.freeThreshold;
clusterThreshold = parameters.clusterThreshold;
meanExtentBirth = ( parameters.priorExtent1 / ( parameters.priorExtent2 - 3 ) ) ^ 2;
measurementsCovariance = parameters.measurementVariance * eye( 2 ) + meanExtentBirth;
minClusterElements = parameters.minClusterElements;

allIndexesNumeric = ( 1:size( measurements, 2 ) )';

freeIndexes = probabilitiesNew >= threshold;
assignedIndexes = probabilitiesNew < threshold;

measurementsFree = measurements( :, freeIndexes );

freeIndexesNumeric = allIndexesNumeric( freeIndexes );
assignedIndexesNumeric = allIndexesNumeric( assignedIndexes );

clusters = getClusters( measurementsFree, measurementsCovariance, clusterThreshold )';

numElements = sum( clusters > 0, 1 );
[ numElements, indexes ] = sort( numElements, 'descend' );
clusters = clusters( :, indexes );

notUsedIndexes = clusters( :, numElements < minClusterElements );
notUsedIndexes = nonzeros( notUsedIndexes( : ) );
notUsedIndexesNumeric = freeIndexesNumeric( notUsedIndexes );
numNotUsed = size( notUsedIndexesNumeric, 1 );

clusters( :, numElements < minClusterElements ) = [  ];
indexesNumericNew = zeros( 0, 1 );
numClusters = size( clusters, 2 );
centralIndexes = zeros( numClusters, 1 );
for cluster = 1:numClusters

    indexes = nonzeros( clusters( :, cluster ) );
    currentMeasurements = measurementsFree( :, indexes );

    currentIndexesNumeric = freeIndexesNumeric( indexes );

    if ( numel( indexes ) > 1 )
        numMeasurements = size( indexes, 1 );
        distanceMatrix = zeros( numMeasurements, numMeasurements );
        for measurement1 = 1:numMeasurements
            for measurement2 = ( measurement1 + 1 ):numMeasurements
                distVector = currentMeasurements( :, measurement1 ) - currentMeasurements( :, measurement2 );

                distanceMatrix( measurement1, measurement2 ) = sqrt( distVector' / measurementsCovariance * distVector );
                distanceMatrix( measurement2, measurement1 ) = distanceMatrix( measurement1, measurement2 );
            end
        end

        distanceVector = sum( distanceMatrix, 2 );
        [ ~, indexes ] = sort( distanceVector, 'descend' );
        currentIndexesNumeric = currentIndexesNumeric( indexes );

    end
    indexesNumericNew = [ currentIndexesNumeric;indexesNumericNew ];

    centralIndexes( 1:cluster ) = centralIndexes( 1:cluster ) + numMeasurements;
end

indexesReordered = [ notUsedIndexesNumeric;indexesNumericNew;assignedIndexesNumeric ];

centralIndexes = centralIndexes + numNotUsed;
centralIndexes = sort( centralIndexes, 'descend' );
end

function [ clusters ] = getClusters( measurements, measurementsCovariance, thresholdProbability )
numMeasurements = size( measurements, 2 );

if ( ~numMeasurements )
    clusters = [  ];
    return ;
end

thresholdDistance = chi2inv( thresholdProbability, 2 );

distanceVector = zeros( ( numMeasurements * ( numMeasurements - 1 ) / 2 + 1 ), 1 );
distanceMatrix = zeros( numMeasurements, numMeasurements );
entry = 1;
for measurement1 = 1:numMeasurements
    for measurement2 = ( measurement1 + 1 ):numMeasurements
        distVector = measurements( :, measurement1 ) - measurements( :, measurement2 );

        entry = entry + 1;
        distanceVector( entry ) = sqrt( distVector' / measurementsCovariance * distVector );

        distanceMatrix( measurement1, measurement2 ) = distanceVector( entry );
        distanceMatrix( measurement2, measurement1 ) = distanceVector( entry );
    end
end

distanceVector = sort( distanceVector );
distanceVector( distanceVector > thresholdDistance ) = [  ];
distance = distanceVector( end  );

clusterNumbers = zeros( numMeasurements, 1 );
clusterId = 1;
for measurement = 1:numMeasurements
    if ( clusterNumbers( measurement ) == 0 )
        clusterNumbers( measurement ) = clusterId;
        clusterNumbers = findNeighbors( measurement, clusterNumbers, clusterId, distanceMatrix, distance );
        clusterId = clusterId + 1;
    end
end
numClusters = clusterId - 1;

maxElements = sum( clusterNumbers == mode( clusterNumbers ) );
clusters = zeros( 0, maxElements );
index = 0;
for cluster = 1:numClusters
    associationTmp = find( clusterNumbers == cluster )';
    numElements = numel( associationTmp );
    if ( numElements <= maxElements )
        index = index + 1;
        clusters( index, : ) = [ zeros( 1, maxElements - numElements ), associationTmp ];
    end
end

end

function [ cellNumbers ] = findNeighbors( index, cellNumbers, cellId, distanceMatrix, distanceThreshold )
numMeasurements = size( distanceMatrix, 2 );

for measurement = 1:numMeasurements
    if ( measurement ~= index && distanceMatrix( measurement, index ) < distanceThreshold && cellNumbers( measurement ) == 0 )
        cellNumbers( measurement ) = cellId;
        cellNumbers = findNeighbors( index, cellNumbers, cellId, distanceMatrix, distanceThreshold );
    end
end

end
