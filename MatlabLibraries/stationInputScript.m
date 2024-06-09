close all
clear
clc
%% setting input structure for the class ssCOD
tic
inputStructure = [];
caseNo = 1;
inputStructure.startTimeVector = [2024, 1, 1, 8, 30, 0.0];
inputStructure.OdDurationInMinutes = 120.0; % in minutes
inputStructure.durationAfterEstimation = 120; % in minutes
inputStructure.plotStepTime = 1.0; % seconds

inputStructure.satelliteMass = 100.0; % in [kg]
inputStructure.orbitAltitude = 700.0; % in [km]
inputStructure.orbitEccentricity = 0.0; % [ndim]
inputStructure.orbitInclination = 89.46; % in [deg]
inputStructure.orbitPregeeArgument = 0.0; % in [deg]
inputStructure.orbitRaan = 180.0; % in [deg]
inputStructure.orbitTrueAnomaly = 0.0; % in [deg]
inputStructure.constellationConfiguration = [420, 15, 1];
inputStructure.underOdSatellite = 367;

inputStructure.underPositioningPoint = [35, 51, 0];
inputStructure.positioningTimeSpanInHours = 0.0167;
inputStructure.positioningIntervals = 1;
inputStructure.maskAngle = 7; % [deg]

% Orbit propagator parameters
inputStructure.propagationMinTimeStep = 0.1;  % s
inputStructure.propagationMaxTimeStep = 3.0;  % stimestep
inputStructure.vectorAbsoluteTolerance = 3.0; % m
inputStructure.vectorRelativeTolerance = 3.0; % m
% other filter parameters
inputStructure.BLSestimatorConvergenceThreshold = 1e-4;
inputStructure.BLSestimatorMaxIterations = 25;
inputStructure.BLSestimatorMaxEvaluations = 35;

inputStructure.showActiveSatellites = false;

inputStructure.stationRangeMaxError = 10.0; % not applicable
inputStructure.stationRangeSigma = 0.01; % not applicable
inputStructure.stationRangeWeight = 0.1; % Not applicable


inputStructure.stationSigma = 1.6; % station measurement sigma
inputStructure.stationBaseWeight = 1.0; % station measurement base weight
inputStructure.stationisTwoWay = true; % TwoWay for those support twoway
inputStructure.stationRefraction = true; % refraction
inputStructure.stationSampleTime = 1.0; % sample time
inputStructure.stationStd = 0.001; % noice sigma

inputStructure.noiseSeed = 10;

inputStructure.gpsSampleTime = 1.0;    % Observation parameter
inputStructure.gpsSigma = 1.8;         % filter parameter
inputStructure.gpsBaseWeight = 1.0;    % filter parameter
inputStructure.gpsPositionStd = 20.34; % GPS Position Error ~ 3-sigma [m]
inputStructure.gpsVelocitystd = 0.2;   % GPS Velocity Error ~ 3-sigma [m/s]

inputStructure.isISL = false;
inputStructure.islLink = 'G01';
inputStructure.islSigma = 1e-4;      % filter parameter
inputStructure.islBaseWeight = 4.0;    % filter parameter
inputStructure.isTwoWay = true;        % Observation parameter
inputStructure.islSampleTime = 1.0;   % Observation parameter
inputStructure.islScheduler = 'Continouos'; % alternative 'InView'
inputStructure.islStd = 1e-4;        % Environment parameter
%%

mainScenario = ssCOD(inputStructure);
mainScenario = mainScenario.feedGroundStations(["Tabriz", "Boushehr", "Bandar Abbas", "Mashhad"]);
underODSatellite = inputStructure.underOdSatellite;
thisConst = inputStructure.constellationConfiguration;
[starConstellation, activeSatellites, selectedOrbitsIndexes, linkList] = ...
    mainScenario.buildISLStar(thisConst(2), thisConst(1)/thisConst(2), ...
    thisConst(3), underODSatellite);
%%
underODOrbit = starConstellation(underODSatellite);
%%
underODPropagator = mainScenario.getPropagator("Get", underODOrbit, "H");
%%
propBuilder = mainScenario.getPropagator("Builder", underODOrbit, "H");
propagatorBuilderArray = javaArray('org.orekit.propagation.conversion.NumericalPropagatorBuilder', 1);
propagatorBuilderArray(1) = propBuilder;
%%
thisEstimator = mainScenario.getPodEstimator(propagatorBuilderArray);
ODStartTime = underODOrbit.getDate();
mainScenario.feedPodMeasurements(thisEstimator, underODOrbit, ODStartTime, mainScenario.duration, "RANGERATE", [])

%%
estimatorOutput = thisEstimator.estimate();
toc
midEstProp = estimatorOutput(1);
midRealProp = underODPropagator;
estEphemerisGenerator = midEstProp.getEphemerisGenerator();
realEphemerisGenerator = midRealProp.getEphemerisGenerator();
midEstProp.propagate(midEstProp.getInitialState().getDate(), mainScenario.endTime);
midRealProp.propagate(mainScenario.startTime, mainScenario.endTime);
estimatedPropagators = estEphemerisGenerator.getGeneratedEphemeris();
realPropagatores = realEphemerisGenerator.getGeneratedEphemeris();
%%
ricposdiff = mainScenario.getRICPositionDifference(realPropagatores, estimatedPropagators);
figureInfo = mainScenario.getFigureTitle(underODSatellite);
%%
thisFigureName = ['thisFigureCaseNo' num2str(caseNo) '.fig'];
thisFigure = mainScenario.getErrorFigure(ricposdiff, figureInfo);
% saveas(thisFigure, thisFigureName)