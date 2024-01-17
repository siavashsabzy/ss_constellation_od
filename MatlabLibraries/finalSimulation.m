close all
clear
clc
%% setting input structure for the class ssCOD
tic
inputStructure = [];
caseNo = 1;
inputStructure.startTimeVector = [2024, 9, 14, 8, 30, 0.0];
inputStructure.OdDurationInMinutes = 15.0; % in minutes
inputStructure.durationAfterEstimation = 90; % in minutes
inputStructure.plotStepTime = 1.0; % seconds

inputStructure.satelliteMass = 100.0; % in [kg]
inputStructure.orbitAltitude = 700.0; % in [km]
inputStructure.orbitEccentricity = 0.0; % [ndim]
inputStructure.orbitInclination = 89.46; % in [deg]
inputStructure.orbitPregeeArgument = 0.0; % in [deg]
inputStructure.orbitRaan = 0.0; % in [deg]
inputStructure.orbitTrueAnomaly = 0.0; % in [deg]
inputStructure.constellationConfiguration = [420, 15, 1];
inputStructure.underOdSatellite = 200;

inputStructure.underPositioningPoint = [35, 51, 0];
inputStructure.positioningTimeSpanInHours = 0.1;
inputStructure.positioningIntervals = 359;
inputStructure.maskAngle = 5; % [deg]

% Orbit propagator parameters
inputStructure.propagationMinTimeStep = 0.1;  % s
inputStructure.propagationMaxTimeStep = 3.0;  % stimestep
inputStructure.vectorAbsoluteTolerance = 3.0; % m
inputStructure.vectorRelativeTolerance = 3.0; % m
% other filter parameters
inputStructure.BLSestimatorConvergenceThreshold = 1e-4;
inputStructure.BLSestimatorMaxIterations = 100;
inputStructure.BLSestimatorMaxEvaluations = 100;

inputStructure.showActiveSatellites = false;

inputStructure.stationRangeMaxError = 0.001;
inputStructure.stationRangeSigma = 0.01;
inputStructure.stationRangeWeight = 0.1;
inputStructure.stationSigma = 0.004; % RangeRate sigma
inputStructure.stationBaseWeight = 1.0;
inputStructure.stationisTwoWay = true;
inputStructure.stationRefraction = true;
inputStructure.stationSampleTime = 1.0;
inputStructure.stationStd = 0.001;

inputStructure.noiseSeed = 1;

inputStructure.gpsSampleTime = 1.0;    % Observation parameter
inputStructure.gpsSigma = 1.8;         % filter parameter
inputStructure.gpsBaseWeight = 1.0;    % filter parameter
inputStructure.gpsPositionStd = 20.34; % GPS Position Error ~ 3-sigma [m]
inputStructure.gpsVelocitystd = 0.2;   % GPS Velocity Error ~ 3-sigma [m/s]

inputStructure.isISL = true;
inputStructure.islLink = 'G01';
inputStructure.islSigma = 1e-4;      % filter parameter
inputStructure.islBaseWeight = 4.0;    % filter parameter
inputStructure.isTwoWay = true;        % Observation parameter
inputStructure.islSampleTime = 0.1;   % Observation parameter
inputStructure.islScheduler = 'Continouos'; % alternative 'InView'
inputStructure.islStd = 1e-4;        % Environment parameter
%%
mainScenario = ssCOD(inputStructure);

underODSatellite = inputStructure.underOdSatellite;
thisConst = inputStructure.constellationConfiguration;

[starConstellation, ~, ~, ~] = ...
    mainScenario.buildISLStar(thisConst(2), thisConst(1)/thisConst(2), ...
    thisConst(3), underODSatellite);

pointLla = inputStructure.underPositioningPoint;
mask = inputStructure.maskAngle;
posSpan = inputStructure.positioningTimeSpanInHours;
posIntervals = inputStructure.positioningIntervals;
% mainScenario.displayLinkedSatellites(pointLla, starConstellation, ...
%     posIntervals, posSpan, mask, 2);

linkedSatelliteStates = mainScenario.getLinkedSatellites(pointLla, ...
    starConstellation, posIntervals, posSpan, mask);

% mainScenario.displayConstellationConfiguration(starConstellation);
% positionErrors = mainScenario.getUserErrorSingleMessage( [31, 51, 0], starConstellation, 14, 0.2, 5.0);
% mainScenario.displayActiveSatellitesConfiguration(activeSatellites, selectedOrbitsIndexes, linkList);
% lgnssEphemeris = mainScenario.getNavigationMessage(starConstellation(1), 900);
% ephemerisError = mainScenario.getLgnssError(lgnssEphemeris, starConstellation(1),...
%     0.25, 1);
%%
allIndexes = [];
for i = 1 : size(linkedSatelliteStates, 2)
    thisStates = linkedSatelliteStates{1,i};
    allIndexes = [allIndexes; thisStates(:,1)];
end
allStates = unique(allIndexes);

step = 1;
NoOfUnderOdOrbits = size(allStates, 1);

clearvars -except mainScenario allStates NoOfUnderOdOrbits j ...
    mainScenarioData inputStructure linkedSatelliteStates
clc

for j = 1 :NoOfUnderOdOrbits
    underODSatellite = allStates(j);
    thisConst = inputStructure.constellationConfiguration;
    [starConstellation, activeSatellites, selectedOrbitsIndexes, linkList] = ...
        mainScenario.buildISLStar(thisConst(2), thisConst(1)/thisConst(2), ...
        thisConst(3), underODSatellite);
    linkedOrbits = [];
    thisLinkList = [underODSatellite];
    underODOrbit = starConstellation(underODSatellite);
    allOrbits = underODOrbit;
    for i = 1:5
        if (linkList(i) ~= -1) && (linkList(i) ~= underODSatellite)
            thisLinkList = [thisLinkList; linkList(i)];
            linkedOrbits = [linkedOrbits;starConstellation(linkList(i))];
            allOrbits = [allOrbits;starConstellation(linkList(i))];
        end
    end
    linkList = []; linkList = thisLinkList';
    linkedPropagators = [];
    underODPropagator = mainScenario.getPropagator("Get", underODOrbit, "H");
    for i = 1:size(linkedOrbits,1)
        linkedPropagators = [linkedPropagators; mainScenario.getPropagator("Get", linkedOrbits(i), "H")];
    end
    propBuilder = [];
    for i  = 1:size(allOrbits,1)
        propBuilder = [propBuilder; mainScenario.getPropagator("Builder", allOrbits(i), "H")];
    end

    thisEstimator = mainScenario.getPodEstimator(propBuilder);
    ODStartTime = underODOrbit(1).getDate();
    mainScenario.feedPodMeasurements(thisEstimator, allOrbits, ODStartTime, mainScenario.duration, "GPS", linkList)
    if mainScenario.isISL
        mainScenario.feedIslMeasurements(thisEstimator, underODPropagator, linkedPropagators, ODStartTime, mainScenario.duration, 'Range', linkList)
        mainScenario.feedIslMeasurements(thisEstimator, underODPropagator, linkedPropagators, ODStartTime, mainScenario.duration, 'Phase', linkList)
    end

    estimatorOutput = thisEstimator.estimate();
    midEstProp = estimatorOutput(1);
    midRealProp = underODPropagator;
    estEphemerisGenerator = midEstProp.getEphemerisGenerator();
    realEphemerisGenerator = midRealProp.getEphemerisGenerator();
    midEstProp.propagate(mainScenario.startTime, mainScenario.endTime);
    midRealProp.propagate(mainScenario.startTime, mainScenario.endTime);
    estimatedPropagator = estEphemerisGenerator.getGeneratedEphemeris();
    realPropagator = realEphemerisGenerator.getGeneratedEphemeris();

    warning('off','MATLAB:rankDeficientMatrix');
    lgnssEphemeris = mainScenario.getNavigationMessageFromEstimatorPropagator(...
        estimatedPropagator, 1800);

    mainScenarioData(j).underODSatellite = underODSatellite;
    mainScenarioData(j).linkList = linkList;
    mainScenarioData(j).underODOrbit = underODOrbit;
    mainScenarioData(j).allOrbits = allOrbits;
    mainScenarioData(j).linkedOrbits = linkedOrbits;
    mainScenarioData(j).linkedPropagators = linkedPropagators;
    mainScenarioData(j).propBuilder = propBuilder;
    mainScenarioData(j).estimator = thisEstimator;
    mainScenarioData(j).estimatorOutput = estimatorOutput;
    mainScenarioData(j).estimatedPropagator = estimatedPropagator;
    mainScenarioData(j).realPropagator = realPropagator;
    mainScenarioData(j).lgnssEphemeris = lgnssEphemeris;
    
    clearvars -except mainScenario allStates NoOfUnderOdOrbits j ...
        mainScenarioData inputStructure linkedSatelliteStates
    clc
end

%%
for i = 1 : size(mainScenarioData, 2)
    mainScenario.saveRICError(mainScenarioData, i, true)
end
%%
positionErrors = mainScenario.getUserPositioningError(...
                inputStructure, linkedSatelliteStates, mainScenarioData);