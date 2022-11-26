% Florian Meyer, 2020

% F. Meyer and J. L. Williams, “Scalable detection and tracking of geometric extended objects,” 
% IEEE Trans. Signal Process., vol. 69, pp. 6283–6298, Oct. 2021.

clear variables; close all; clc; addpath './_common'; rng(1);


% parameters of simulated scenario
numSteps = 50;
numTargets = 5;
meanTargetDimension = 3;
startRadius = 75;
startVelocity = 10;


% main parameters of the statistical model
parameters.scanTime = .2;
parameters.accelerationDeviation = 1;
parameters.survivalProbability = 0.99;
parameters.meanBirths = .01;
parameters.surveillanceRegion = [[-200; 200] [-200; 200]];
parameters.measurementVariance = 1^2;
parameters.meanMeasurements = 8;
parameters.meanClutter = 10;


% prior distribution parameters
parameters.priorVelocityCovariance = diag([10^2;10^2]);
parameters.priorExtent2 = 100;                       
parameters.priorExtent1 = [[meanTargetDimension 0];[0 meanTargetDimension]]*(parameters.priorExtent2-3);
parameters.degreeFreedomPrediction = 20000;

% censoring and measurement reordering parameters
parameters.freeThreshold = 0.9;
parameters.clusterThreshold = 0.9;
parameters.minClusterElements = 1;


% sampling parameters
parameters.numParticles = 5000;
parameters.regularizationDeviation = 0;


% detection and pruning parameters
parameters.detectionThreshold = .5;
parameters.thresholdPruning = 10^(-3);
parameters.minimumTrackLength = 1;


% message passing parameters
parameters.numOuterIterations = 2;



% generate true start states
[startStates,startMatrixes] = getStartStates(numTargets,startRadius,startVelocity,parameters);
appearanceFromTo = [[3;83],[3;83],[6;86],[6;86],[9;89],[9;89],[12;92],[12;92],[15;95],[15;95]];


% generate true track
[targetTracks,targetExtents] = generateTracksUnknown(parameters,startStates,startMatrixes,appearanceFromTo,numSteps);


% generate measurements
measurements = generateClutteredMeasurements(targetTracks, targetExtents, parameters);


% perform graph-based extended object tracking (EOT)
tic
[estimatedTracks,estimatedExtents] = eotEllipticalShape(measurements,parameters);
toc


% show results
mode = 1; %hit ``space'' to start visualization; set mode=0 for final result and mode=2 to frame-by-frame.
showResults(targetTracks,targetExtents,estimatedTracks,estimatedExtents,measurements,[-150 150 -150 150],mode);
