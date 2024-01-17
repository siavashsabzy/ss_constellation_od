close all
clear
clc
caseObjects = loadCaseStudies();
for caseNo = 1:size(caseObjects,1)
    mainScenario = caseObjects(caseNo);
    underODSatellite = mainScenario.thisSatellite(4);
    numPlane = mainScenario.thisSatellite(2);
    numSats = mainScenario.thisSatellite(1)/mainScenario.thisSatellite(2);
    phase = mainScenario.thisSatellite(3);
    [starConstellation, activeSatellites, selectedOrbitsIndexes, linkList] = ...
        mainScenario.buildISLStar(numPlane, numSats, phase, underODSatellite);

    % mainScenario.displayConstellationConfiguration(starConstellation);
    % mainScenario.displayActiveSatellitesConfiguration(activeSatellites, selectedOrbitsIndexes, linkList);

    %%
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
    thisFigureName = ['thisFigureCaseNo' num2str(caseNo) '.fig'];
    thisFigure = mainScenario.getErrorFigure(ricposdiff, figureInfo);
    saveas(thisFigure, thisFigureName)
    pause(1.0)
    close all
end