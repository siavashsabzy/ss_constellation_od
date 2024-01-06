close all
clear
clc
warning('off','MATLAB:rankDeficientMatrix');
%%
tic
try
    mainScenario = ssCOD();
catch
    % ssCOD.addOrekitLibraries();
    mainScenario = ssCOD();
end
%mainScenario.plotReferenceOrbitTrajectory(5,60);
underODSatellite = 10;
[starConstellation, activeSatellites, selectedOrbitsIndexes, linkList] = ...
    mainScenario.buildISLStar(11, 23, 1, underODSatellite);
% positionErrors = mainScenario.getUserErrorSingleMessage( [31, 51, 0], starConstellation, 14, 0.2, 5.0);
% [GDOP, PDOP, TDOP, HDOP, VDOP] = mainScenario.getGlobeDops(100000, starConstellation, 5, true);
[GDOP, PDOP, TDOP, HDOP, VDOP] = mainScenario.getGlobeDopsHistory(1000, starConstellation, 5,...
    mainScenario.referenceOrbit.getKeplerianPeriod()/3600, 10);
% mainScenario.displayConstellationConfiguration(starConstellation);
% mainScenario.displayActiveSatellitesConfiguration(activeSatellites, selectedOrbitsIndexes, linkList);
% lgnssEphemeris = mainScenario.getNavigationMessage(starConstellation(1), 900);
% ephemerisError = mainScenario.getLgnssError(lgnssEphemeris, starConstellation(1),...
%     0.25, 1);
%%
% ecef = mainScenario.eph2xyz(lgnssEphemeris, 0);
toc
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
%%
linkedPropagators = [];
underODPropagator = mainScenario.getPropagator("Get", underODOrbit, "H");
for i = 1:size(linkedOrbits,1)
    linkedPropagators = [linkedPropagators; mainScenario.getPropagator("Get", linkedOrbits(i), "H")];
end
%%
propBuilder = [];
for i  = 1:size(allOrbits,1)
    propBuilder = [propBuilder; mainScenario.getPropagator("Builder", allOrbits(i), "H")];
end
%%
thisEstimator = mainScenario.getPodEstimator(propBuilder);
ODStartTime = underODOrbit(1).getDate();
if mainScenario.isISL
    mainScenario.feedIslMeasurements(thisEstimator, underODPropagator, linkedPropagators, ODStartTime, mainScenario.duration, 'Range', linkList)
    mainScenario.feedIslMeasurements(thisEstimator, underODPropagator, linkedPropagators, ODStartTime, mainScenario.duration, 'Phase', linkList)
end
mainScenario.feedPodMeasurements(thisEstimator,allOrbits, ODStartTime, mainScenario.duration, "GPS", linkList)

%%
estimatorOutput = thisEstimator.estimate();
estimatedOrbits = [];
estimatedPropagators = [];
realPropagatores = [];
for indx = 1:size(propBuilder,1)
    thisOrbit = estimatorOutput(indx).getInitialState().getOrbit();
    estimatedOrbits = [estimatedOrbits; thisOrbit];
    midEstProp = mainScenario.getPropagator("Get", thisOrbit, "H");
    midRealProp = mainScenario.getPropagator("Get", allOrbits(indx), "H");
    estEphemerisGenerator = midEstProp.getEphemerisGenerator();
    realEphemerisGenerator = midRealProp.getEphemerisGenerator();
    midEstProp.propagate(mainScenario.startTime, mainScenario.endTime);
    midRealProp.propagate(mainScenario.startTime, mainScenario.endTime);
    estimatedPropagators = [estimatedPropagators; estEphemerisGenerator.getGeneratedEphemeris()];
    realPropagatores = [realPropagatores; realEphemerisGenerator.getGeneratedEphemeris()];
    clear thisOrbit
end
%%
Sat = 1;
ricposdiff = mainScenario.getRICPositionDifference(realPropagatores(Sat), estimatedPropagators(Sat));
figureInfo = mainScenario.getFigureTitle(thisLinkList(Sat));
%%
mainScenario.getErrorFigure(ricposdiff, figureInfo)