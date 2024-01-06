import orekit

orekit.initVM()
from orekit import JArray, JArray_double
from org.orekit.models.earth import ReferenceEllipsoid
from org.orekit.attitudes import NadirPointing
from org.orekit.propagation.events import ElevationDetector, InterSatDirectViewDetector
from org.orekit.propagation.events.handlers import ContinueOnEvent
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.gravity import ThirdBodyAttraction
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient, IsotropicRadiationCNES95Convention
from org.orekit.forces.radiation import SolarRadiationPressure
from org.orekit.forces.gravity import Relativity, OceanTides, SolidTides
from org.orekit.models.earth.atmosphere.data import MarshallSolarActivityFutureEstimation
from org.orekit.models.earth.atmosphere import DTM2000
from org.orekit.forces.drag import IsotropicDrag
from org.orekit.forces.drag import DragForce
from org.orekit.estimation.leastsquares import BatchLSEstimator, PythonBatchLSObserver
from org.hipparchus.optim.nonlinear.vector.leastsquares import GaussNewtonOptimizer
from org.hipparchus.linear import QRDecomposer
from org.orekit.propagation.conversion import NumericalPropagatorBuilder
from org.orekit.estimation.measurements import GroundStation, Range, RangeRate, ObservableSatellite, AngularAzEl, AngularRaDec, \
    Position, PV, InterSatellitesRange
from org.orekit.estimation.measurements.gnss import OneWayGNSSRange, InterSatellitesPhase
from org.orekit.estimation.measurements.generation import Generator, RangeBuilder, RangeRateBuilder, AngularAzElBuilder, \
    EventBasedScheduler, SignSemantic, AngularRaDecBuilder, InterSatellitesRangeBuilder, InterSatellitesPhaseBuilder, \
    ContinuousScheduler, OneWayGNSSRangeBuilder, GatheringSubscriber, PVBuilder
       

from org.orekit.gnss import  Frequency
# from org.orekit.estimation.measurements.generation import RangeBuilder, RangeRateBuilder
from org.orekit.frames import FramesFactory, TopocentricFrame, ITRFVersion, LOFType
from org.orekit.time import AbsoluteDate, TimeScalesFactory, FixedStepSelector, BurstSelector
from org.orekit.utils import Constants, IERSConventions, PVCoordinates, TimeStampedPVCoordinates
from org.orekit.bodies import GeodeticPoint, CelestialBodyFactory
from org.orekit.orbits import KeplerianOrbit, CartesianOrbit, PositionAngleType
from org.orekit.propagation import SpacecraftState
# from org.orekit.propagation.analytical import KeplerianPropagator
# from org.orekit.propagation.sampling import OrekitFixedStepHandler
from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation.numerical import NumericalPropagator
# from org.hipparchus.random import CorrelatedRandomVectorGenerator, NormalizedRandomGenerator
# from org.hipparchus.linear import RealMatrix
from org.hipparchus.geometry.euclidean.threed import Vector3D, SphericalCoordinates
from org.hipparchus.random import Well19937a, GaussianRandomGenerator, CorrelatedRandomVectorGenerator
from org.hipparchus.linear import MatrixUtils, Array2DRowRealMatrix

# Atmospheric Models
from org.orekit.models.earth import EarthITU453AtmosphereRefraction
# Tropospheric Models
from org.orekit.models.earth.troposphere import SaastamoinenModel

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
from org.orekit.propagation.analytical.gnss import GNSSPropagator, GNSSPropagatorBuilder
from org.orekit.propagation.analytical.gnss.data import CommonGnssData, GNSSOrbitalElements, AbstractNavigationMessage
from org.orekit.propagation.conversion import PropagatorConverter, FiniteDifferencePropagatorConverter, AbstractPropagatorConverter
from org.orekit.estimation.iod import IodGooding



setup_orekit_curdir()
import math
import numpy as np
import random
import plotly.express as px
import plotly.io as pio
import plotly.graph_objs as go
import pandas as pd
from bs4 import BeautifulSoup
from PIL import Image
from datetime import datetime
import webbrowser

webbrowser.register('firefox', None, webbrowser.BackgroundBrowser("C:\\Program Files\\Mozilla Firefox\\firefox.exe"))
pio.renderers.default = "firefox"

class ssEnvironmentBuilder():

    def __init__(self):
        super().__init__()       
        self.utc = TimeScalesFactory.getUTC()
        self.startTime = AbsoluteDate(2020, 8, 13, 19, 0, 0.0, self.utc)
        self.duration = 10 / 60
        self.durationAfterEstimation = 45.0 / 60.0 
        self.endEstimation = self.startTime.shiftedBy(self.duration * 3600.0)
        self.endTime = self.endEstimation.shiftedBy(self.durationAfterEstimation * 3600.0)
        self.plotStepTime = 1.0
        #self.endTime = AbsoluteDate(2020, 8, 14, 19, 0, 0.0, self.utc)
        self.RAD2DEG = (180 / np.pi)
        self.DEG2RAD = (np.pi / 180)
        self.nominalAltitudes = [600, 800] # range

        self.showActiveSatellites = False


        self.rangeMaxError = 0.001

        self.rangeSigma = 0.01
        self.rangeWeight = 0.1
        self.stationSigma=0.00400 # RangeRate sigma
        self.stationBaseWeight=1.0
        self.stationisTwoWay = True
        self.stationRefractio = True
        self.stationSampleTime = 1.0
        self.stationStd = 0.001

        self.noiseSeed = 1


        self.isISL = False
        self.islSigma = 0.005       # filter parameter
        self.islBaseWeight = 4.0    # filter parameter
        self.isTwoWay = True        # Observation parameter
        self.islSampleTime = 0.10   # Observation parameter
        self.islScheduler = 'Continuous' # alternative 'InView'
        self.islStd = 0.0001        # Environment parameter

        self.gpsSampleTime = 1.0    # Observation parameter
        self.gpsSigma = 1.8         # filter parameter
        self.gpsBaseWeight = 1.0    # filter parameter
        self.gpsPositionStd = float(100) # GPS Position Error ~ 3-sigma: 100 [m]
        self.gpsVelocitystd = self.gpsPositionStd / 20 # GPS Velocity Error
        self.GPSMaxBurstSize = 1    # maximum number of selected dates in a burst
        self.GPSHighRateStep = 1.0 * self.gpsSampleTime # step between two consecutive dates within a burst (s)
        self.GPSBurstPeriod  = 1.0 * self.gpsSampleTime # period between the start of each burst (s)




        self.groundStationsList = []
        self.numberOfGroundStation = 0
        self.spacecraftsList = []
        self.propagatorsList = []
        self.referenceOrbit = []
        
        self.eciFrame = FramesFactory.getGCRF()
        self.ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
        self.wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(self.ecefFrame)
        self.moon = CelestialBodyFactory.getMoon()
        self.sun = CelestialBodyFactory.getSun()  
        # Orbit propagator parameters
        self.propagationMinTimeStep = 0.5  # s
        self.propagationMaxTimeStep = 2.0  # stimestep
        self.vectorAbsoluteTolerance = 3.0  # m
        self.vectorRelativeTolerance = 3.0 # m
        # other filter parameters
        self.BLSestimatorConvergenceThreshold = 1e-5
        self.BLSestimatorMaxIterations = 25
        self.BLSestimatorMaxEvaluations = 35



        self.elevationMask = 10.0 * np.pi / 180.0 # radian(s)

        self.earthMapColorScale =[[0.0, 'rgb(30, 59, 117)'],
                 [0.1, 'rgb(46, 68, 21)'],
                 [0.2, 'rgb(74, 96, 28)'],
                 [0.3, 'rgb(115,141,90)'],
                 [0.4, 'rgb(122, 126, 75)'],
                 [0.6, 'rgb(122, 126, 75)'],
                 [0.7, 'rgb(141,115,96)'],
                 [0.8, 'rgb(223, 197, 170)'],
                 [0.9, 'rgb(237,214,183)'],
                 [1.0, 'rgb(255, 255, 255)']]
        

        # referenceDate = AbsoluteDate.J2000_EPOCH  # reference date
        # mjdUtcEpoch = AbsoluteDate(1858, 11, 17, 0, 0, 0.0, utc)

    def feedGroundStations(self, names):
        '''
        this method is designed to define a collection of ground stations to the 
        environment model. The method only take one or a list of following stations:
        Boushehr, Tehran, Aleshtar, Bandar Abbas, Mashhad, Tabriz, Sistan
        
        Each time you run the method, the collection will be intialized.
        '''

        addedGroundStations = []
        for i in range(len(names)):
            if names[i] == "Boushehr":
                self.addGroundStation(30.0, 48.5, 0.0, "Boushehr")
                addedGroundStations.append(names[i])
            elif names[i] == "Tehran":
                self.addGroundStation(35.0, 51.0, 0.0, "Tehran")
                addedGroundStations.append(names[i])
            elif names[i] == "Aleshtar":
                self.addGroundStation(33.0, 48.0, 1500.0, "Aleshtar")
                addedGroundStations.append(names[i])
            elif names[i] == "Bandar Abbas":
                self.addGroundStation(25.0, 59.0, 0.0, "Bandar Abbas")
                addedGroundStations.append(names[i])
            elif names[i] == "Mashhad":
                self.addGroundStation(35.0, 58.0, 0.0, "Mashhad")
                addedGroundStations.append(names[i])
            elif names[i] == "Tabriz":
                self.addGroundStation(38.0, 46.0, 0.0, "Tabriz")
                addedGroundStations.append(names[i])
            elif names[i] == "Sistan":
                self.addGroundStation(27.0, 63.0, 0.0, "Sistan")
                addedGroundStations.append(names[i])
            else:
                pass
        print(" ")
        print("added Ground Stations ({}): {}".format(len(addedGroundStations), addedGroundStations))
            
            
            
            

    def addReferenceOrbit(self, altitude, eccentricity, inclination, raan, aop, trueAnomaly,
                           anomalyType=PositionAngleType.TRUE):
        '''
        this method is designed to define a reference orbit to the environment, Based on this orbit the 
        following constellation will be built.
        
        Each time you run the method, the orbit will be intialized.
        '''

        self.referenceOrbit = []
        self.referenceOrbit = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + altitude,  # Semimajor Axis (m)
                                      eccentricity,  # Eccentricity
                                      inclination * np.pi / 180,  # Inclination (radians)
                                      aop * np.pi / 180,  # Perigee argument (radians)
                                      raan * np.pi / 180,  # Right ascension of ascending node (radians)
                                      trueAnomaly * np.pi / 180,  # Anomaly (rad/s)
                                      anomalyType,  # Sets which type of anomaly we use
                                      self.eciFrame, # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                                      self.startTime,  # Sets the date of the orbital parameters
                                      Constants.WGS84_EARTH_MU)  # Sets the central attraction coefficient (m³/s²)
        
    
    def performGooding(self, radecs, ranges):
        '''
        this method is designed to select proper measurements and perform an iod as initial guesses for pod
            as inputs this function takes lists of rightAscentions-declination, ranges and validated altitude range
        if successful this method gives a keplerian orbit
        '''
        # select measurements
        IodFunction = IodGooding(Constants.WGS84_EARTH_MU)
        orbit = []
        lowerBand = self.nominalAltitudes[0] * 1000.0 + Constants.WGS84_EARTH_EQUATORIAL_RADIUS
        upperBand = self.nominalAltitudes[1] * 1000.0 + Constants.WGS84_EARTH_EQUATORIAL_RADIUS
        try:
            initialRADEC = radecs[0]
            initialRange = ranges[0]
            lastRADEC = radecs[-1]
            lastRange = ranges[-1]
            middleRADEC = radecs[int(len(radecs)/2)]

            if (  (initialRADEC.getDate().compareTo(initialRange.getDate()) > 0.01)    or 
                    (lastRADEC.getDate().compareTo(lastRange.getDate()) > 0.01)     ):
                pass
            else:
                orbitGuess = IodFunction.estimate(self.eciFrame, initialRADEC, middleRADEC, lastRADEC,
                                    initialRange.getObservedValue()[0], lastRange.getObservedValue()[0], 0, True)
                if (orbitGuess.getA() < lowerBand) or (orbitGuess.getA() > upperBand):
                    pass
                else:
                    orbit = orbitGuess
        except:
            pass    
        return orbit
    
    def performIod(self, propagator, groundStationsList, startTime , duration=1.5):
        orbit = []
        while orbit == []:
            for i in range(len(groundStationsList)):
                radecs = []
                ranges = []
                station = groundStationsList[i]
                stationName = groundStationsList[i].getBaseFrame().getName()   
                altitude = groundStationsList[i].getBaseFrame().getPoint().getAltitude()
                radecs = self.generateMeasurements(propagator, station, stationName, altitude,
                                        self.elevationMask, "RADEC",
                                        startTime, duration, self.eciFrame)
                ranges = self.generateMeasurements(propagator, station, stationName , altitude,
                                        self.elevationMask, "RANGE",
                                        startTime, duration, self.eciFrame)
                
                orbit = self.performGooding(radecs, ranges)

                if orbit != []:
                    break
                else:
                    orbit = []
            startTime = startTime.shiftedBy(duration*3600.0)

        return orbit


    def getPodEstimator(self, propagatorBuilder):
        matrixDecomposer = QRDecomposer(1e-11)
        optimizer = GaussNewtonOptimizer(matrixDecomposer, False)
        satEstimator = BatchLSEstimator(optimizer, propagatorBuilder)
        satEstimator.setParametersConvergenceThreshold(self.BLSestimatorConvergenceThreshold)
        satEstimator.setMaxIterations(self.BLSestimatorMaxIterations)
        satEstimator.setMaxEvaluations(self.BLSestimatorMaxEvaluations)
        satEstimator.setObserver(MyObserver())
        return satEstimator
    
    def feedPodMeasurements(self, estimator, allOrbits, startTime, duration, measurementsType, linkList):
        
        if measurementsType == "GPS":
            meas = self.generateMeasurements(allOrbits, [], [], 0.0,
                                        self.elevationMask, "GPS",
                                            startTime, duration, self.eciFrame, linkList)
            i = 0
            for m in meas:
                estimator.addMeasurement(m)
                i += 1
            print("{} measurements added to estimator".format(i))
        else:
            for i in range(self.numberOfGroundStation):
                station = self.groundStationsList[i]
                stationName = self.groundStationsList[i].getBaseFrame().getName()   
                altitude = self.groundStationsList[i].getBaseFrame().getPoint().getAltitude()
                meas = self.generateMeasurements(allOrbits, station, stationName, altitude,
                                            self.elevationMask, "RANGERATE",
                                                startTime, duration, self.eciFrame, linkList)
                i = 0
                for m in meas:
                    estimator.addMeasurement(m)
                    i += 1
                print("{} measurements added to estimator, for {} station".format(i, stationName))


    def feedIslMeasurements(self, estimator, underODPropagator, linkedpropagators, startTime, duration, islType, linkList):
        initIslTime = startTime
        finalIslTime = initIslTime.shiftedBy(duration * 3600.0)
        mn = 0
        fixed_step_selector = FixedStepSelector(self.islSampleTime , self.utc)
        small = 1e-10
        seed = self.noiseSeed
        random_generator = Well19937a(int(seed))
        gaussian_generator = GaussianRandomGenerator(random_generator)
        covariance = MatrixUtils.createRealDiagonalMatrix([self.islStd**2, self.islStd**2])
        noise_source = CorrelatedRandomVectorGenerator(covariance, float(small), gaussian_generator)
        gatherer = GatheringSubscriber()
        generator = Generator()
        generator.addSubscriber(gatherer)
        local = underODPropagator
        localEphemerisGenerator = local.getEphemerisGenerator()
        local.propagate(initIslTime , finalIslTime)
        localBoundedPropagator = localEphemerisGenerator.getGeneratedEphemeris()
        satellite = generator.addPropagator(localBoundedPropagator)
        Rsatellites = []
        for i in range(len(linkedpropagators)):
            remote = linkedpropagators[i]
            remoteEphemerisGenerator = remote.getEphemerisGenerator()
            remote.propagate(initIslTime , finalIslTime)
            remoteBoundedPropagator = remoteEphemerisGenerator.getGeneratedEphemeris()
            Rsatellites.append(generator.addPropagator(remoteBoundedPropagator))
        for j in range(len(Rsatellites)):
            if islType == 'Range':
                Builder = InterSatellitesRangeBuilder(noise_source, satellite, 
                                                                    Rsatellites[j], self.isTwoWay, 
                                                                    self.islSigma, self.islBaseWeight)
            else:
                Builder = InterSatellitesPhaseBuilder(noise_source, satellite, 
                                                                Rsatellites[j], Frequency.G01.getWavelength(), 
                                                                self.islSigma, self.islBaseWeight)
            if self.islScheduler == 'InView':
                interDetector = InterSatDirectViewDetector( self.wgs84Ellipsoid , remoteBoundedPropagator)
                scheduler = EventBasedScheduler(Builder, fixed_step_selector, underODPropagator,
                                             interDetector,
                                            SignSemantic.FEASIBLE_MEASUREMENT_WHEN_POSITIVE)
                generator.addScheduler(scheduler)
            else:
                generator.addScheduler(ContinuousScheduler(Builder, fixed_step_selector))
            generator.generate(initIslTime, finalIslTime)
            measurements = gatherer.getGeneratedMeasurements()
            for measObject in measurements:
                mn += 1
                if islType == 'Range':
                    meas = InterSatellitesRange.cast_(measObject)
                else:
                    meas = InterSatellitesPhase.cast_(measObject)
                estimator.addMeasurement(meas)
            print("{} measurements feeded to satellite #{} estimator from ISL {} observations of satellite #{}".format(mn,
                   linkList[0], islType, linkList[j+1]))
         
    def addGroundStation(self, latitude, longitude, altitude, name):
        stationPoint = GeodeticPoint(latitude * self.DEG2RAD, longitude * self.DEG2RAD, altitude)
        stationFrame = TopocentricFrame(self.wgs84Ellipsoid, stationPoint, name)
        groundStation = GroundStation(stationFrame)
        self.groundStationsList.append(groundStation)
        self.numberOfGroundStation += 1

    def addSpacecraft(self, spacecraft):
        self.spacecraftsList.append(spacecraft)

    def displayActiveSatellitesConfiguration(self, activeSatellites, selectedOrbitsIndexes, linkList):
        earhtMarbleFigure = self.getBlueMarbleFigure('earth.jpeg')
        epochECI2ECEF = self.eciFrame.getTransformTo(self.ecefFrame, self.startTime)
        for i in range(len(selectedOrbitsIndexes)):
            if selectedOrbitsIndexes[i] == linkList[0, 0]:
                mainSat = i
        mainSatEci = activeSatellites[mainSat].getPVCoordinates().getPosition()
        mainSatECEF = epochECI2ECEF.getRotation().applyTo(mainSatEci)

        constellationX = []
        constellationY = []
        constellationZ = []
        for j in range(len(activeSatellites)):
            if j == mainSat:
                pass
            else:
                eciPos = activeSatellites[j].getPVCoordinates().getPosition()
                truePositionECEF = epochECI2ECEF.getRotation().applyTo(eciPos)
                constellationX.append(truePositionECEF.getX())
                constellationY.append(truePositionECEF.getY())
                constellationZ.append(truePositionECEF.getZ())
        earhtMarbleFigure.add_trace(go.Scatter3d(x = constellationX, 
                                            y = constellationY, 
                                            z = constellationZ,  mode='markers'))
        earhtMarbleFigure.add_trace(go.Scatter3d(x = [mainSatECEF.getX()], 
                                            y = [mainSatECEF.getY()], 
                                            z = [mainSatECEF.getZ()],  mode='markers'))
        earhtMarbleFigure.update_coloraxes(showscale=False)
        earhtMarbleFigure.update_layout(showlegend=False)
        pio.show(earhtMarbleFigure)

    def getInertialPositionDifference(self, prop1, prop2, currentDateTime, endTime, step):
        positionResidual = pd.DataFrame(columns=['X', 'Y', 'Z','Norm'])
        while currentDateTime.compareTo(endTime) <= 0:
            truePosition = prop1.propagate(currentDateTime).getPVCoordinates().getPosition()
            estimatedPosition = prop2.propagate(currentDateTime).getPVCoordinates().getPosition()
            positionDifference = truePosition.subtract(estimatedPosition)
            positionResidual.loc[absolutedate_to_datetime(currentDateTime)] = \
                [positionDifference.getX(), positionDifference.getY(), positionDifference.getZ(), positionDifference.getNorm()]
            currentDateTime = currentDateTime.shiftedBy(step)
        return positionResidual

    def getRICPositionDifference(self, prop1, prop2, currentDateTime, endTime, step):
        positionResidual = pd.DataFrame(columns=['X', 'Y', 'Z','Norm'])
        while currentDateTime.compareTo(endTime) <= 0:
            truePosition = prop1.propagate(currentDateTime).getPVCoordinates().getPosition()
            estimatedPosition = prop2.propagate(currentDateTime).getPVCoordinates().getPosition()
            positionDifference = truePosition.subtract(estimatedPosition)
            errorInLVLH = LOFType.LVLH.rotationFromInertial(prop1.propagate(currentDateTime).getPVCoordinates()).applyTo(positionDifference)
            positionResidual.loc[absolutedate_to_datetime(currentDateTime)] = \
                [errorInLVLH.getX(), errorInLVLH.getY(), errorInLVLH.getZ(), errorInLVLH.getNorm()]
            currentDateTime = currentDateTime.shiftedBy(step)
        return positionResidual
    
    def getResidualFigure(self, dataFrame, additionalInfoStrings):
        trace = go.Scattergl(
                    x=dataFrame.index, y=dataFrame['Norm'],
                    mode='markers',
                    name='PositionNormError'
                )
        layout = go.Layout(title='Position residuals: ' + additionalInfoStrings, 
                   xaxis=dict(title='Datetime UTC'), 
                   yaxis=dict(title='Position residual (m)') )

        fig = dict(data=trace, layout=layout)
        return fig

    def generateMeasurements(self, allOrbits, station, station_name, altitude, elevAngle, meas_type, t0, duration, eciFrame, linkList, mass = 100.0):
        two_way = self.stationisTwoWay
        withRefraction= self.stationRefractio
        step = self.stationSampleTime
        seed = self.noiseSeed
        tf = t0.shiftedBy(3600.0 * duration)
        measurementslist = []
        for proIndex in range(len(allOrbits)):
            gatherer = GatheringSubscriber()
            generator = Generator()
            generator.addSubscriber(gatherer)
            for proIndex2 in range(len(allOrbits)): 
                if proIndex2 == proIndex:
                    satellite = generator.addPropagator(mainScenario.getPropagator("Get", allOrbits[proIndex2], mass, "H"))
                else:
                    generator.addPropagator(mainScenario.getPropagator("Get", allOrbits[proIndex2], mass, "H"))
            fixed_step_selector = FixedStepSelector(step, TimeScalesFactory.getUTC())
            sigma = self.stationSigma
            base_weight = self.stationBaseWeight
            small = 1e-10
            random_generator = Well19937a(int(seed))
            gaussian_generator = GaussianRandomGenerator(random_generator)
            covariance = MatrixUtils.createRealDiagonalMatrix([float(sigma * sigma), float(sigma * sigma)])
            noise_source = CorrelatedRandomVectorGenerator(covariance, float(small), gaussian_generator)
            if meas_type == 'RANGE':
                builder = RangeBuilder(noise_source, station, two_way, sigma, base_weight, satellite)
            elif meas_type == 'RANGERATE':
                builder = RangeRateBuilder(noise_source, station, two_way, sigma, base_weight, satellite)
            elif meas_type == 'AZEL':
                sigmaAzEl = JArray_double([sigma, sigma])
                weightAzEl = JArray_double([base_weight, base_weight])
                builder = AngularAzElBuilder(noise_source, station, sigmaAzEl, weightAzEl, satellite)
            elif meas_type == 'RADEC':
                sigmaRaDec = JArray_double([sigma, sigma])
                weightRaDec = JArray_double([base_weight, base_weight])
                builder = AngularRaDecBuilder(noise_source, station, eciFrame, sigmaRaDec, weightRaDec, satellite)
            elif meas_type == 'GPS':
                sigmaP = self.gpsPositionStd * 2
                sigmaV = self.gpsVelocitystd * 2
                GPSCov = MatrixUtils.createRealDiagonalMatrix([sigmaP, sigmaP, sigmaP, sigmaV, sigmaV, sigmaV])
                GPSNoise = CorrelatedRandomVectorGenerator(GPSCov, float(small), gaussian_generator)
                GPSBuilder = PVBuilder(GPSNoise, self.gpsSigma, self.gpsSigma*0.05, self.gpsBaseWeight, satellite)
                GPSsample = BurstSelector(self.GPSMaxBurstSize, self.GPSHighRateStep, self.GPSBurstPeriod, TimeScalesFactory.getUTC())
                GPSScheduler = ContinuousScheduler(GPSBuilder, GPSsample)
            if meas_type == 'GPS':
                generator.addScheduler(GPSScheduler)
            else:
                if withRefraction:
                    elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(
                        elevAngle).withRefraction(EarthITU453AtmosphereRefraction(altitude)).withHandler(ContinueOnEvent())
                else:
                    elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(elevAngle).withHandler(
                        ContinueOnEvent())
                scheduler = EventBasedScheduler(builder, fixed_step_selector, 
                                                mainScenario.getPropagator("Get", allOrbits[proIndex2], mass, "H"), 
                                                elevation_detector,
                                                SignSemantic.FEASIBLE_MEASUREMENT_WHEN_POSITIVE)
                generator.addScheduler(scheduler)
            generator.generate(t0, tf)
            measurements = gatherer.getGeneratedMeasurements()
            j = 0
            measTime = []
            for measObject in measurements:
                j += 1
                if meas_type == 'RANGE':
                    meas = Range.cast_(measObject)
                elif meas_type == 'RANGERATE':
                    meas = RangeRate.cast_(measObject)
                elif meas_type == 'AZEL':
                    meas = AngularAzEl.cast_(measObject)
                elif meas_type == 'RADEC':
                    meas = AngularRaDec.cast_(measObject)
                elif meas_type == 'GPS':
                    meas = PV.cast_(measObject)
                else:
                    raise ValueError("unrecognized measurement type")
                measurementslist.append(meas)
                if j == 1:
                    measTime.append(meas)
                if two_way:
                    way_str = "TWOWAY"
                else:
                    way_str = "ONEWAY"
            try:
                if meas_type == 'GPS':
                    measTime.append(meas)
                    print("{} number of GPS observations are generated for satellite #{}, first at: {} and last at: {}".format(j, 
                            linkList[proIndex], measTime[0].getDate().toString(), measTime[1].getDate().toString()))
                else:
                    measTime.append(meas)
                    print("{} number of {} observations are generated for satellite #{}, first at: {} and last at: {}, - station {}".format(j, 
                        meas_type, linkList[proIndex], measTime[0].getDate().toString(), measTime[1].getDate().toString(), station_name))
            except:
                pass
        return measurementslist

    def getRangeMeasurements(self, propagator, absoluteDateTime, elevationMask):
        currentInertialCoordinates = propagator.propagate(absoluteDateTime).getPVCoordinates()
        rangeMeasurementsList = []
        for groundStation in self.groundStationsList:
            elevation = groundStation.getBaseFrame().getElevation(
                currentInertialCoordinates.getPosition(), propagator.getFrame(), absoluteDateTime)
            if ((elevation * self.RAD2DEG) > elevationMask):
                realRange = groundStation.getBaseFrame().getRange(currentInertialCoordinates.getPosition(),
                                                                  propagator.getFrame(), absoluteDateTime)
                MeasuredRange = random.uniform(-self.rangeMaxError, self.rangeMaxError) + realRange

                Measurement = Range(groundStation, self.isTwoWay, absoluteDateTime, MeasuredRange,
                                    self.rangeSigma, self.rangeWeight, ObservableSatellite(0))

                rangeMeasurementsList.append(Measurement)

        return rangeMeasurementsList

    def getPropagator(self, propagatorType, initialOrbit, mass, level):

        eciFrame = self.eciFrame
        ecefFrame = self.ecefFrame
        minStep = self.propagationMaxTimeStep
        maxStep = self.propagationMaxTimeStep
        vecAbsoluteTolerance = self.vectorAbsoluteTolerance
        vecRelativeTolerance = self.vectorRelativeTolerance
        moon = self.moon
        sun = self.sun
        wgs84Ellipsoid = self.wgs84Ellipsoid
        nadirPointing = NadirPointing(eciFrame, wgs84Ellipsoid)

        initialCartesianOrbit = CartesianOrbit(SpacecraftState(initialOrbit, mass).getPVCoordinates(eciFrame),
                                               eciFrame, wgs84Ellipsoid.getGM())

        if (propagatorType == 'Builder'):
            integratorBuilder = DormandPrince853IntegratorBuilder(minStep, maxStep, vecAbsoluteTolerance)

            satPropagator = NumericalPropagatorBuilder(initialCartesianOrbit, integratorBuilder, PositionAngleType.TRUE,
                                                       1.0)
            satPropagator.setAttitudeProvider(nadirPointing)
            satPropagator.setMass(mass)
        else:
            thisintegrator = DormandPrince853Integrator(minStep,
                                                        maxStep,
                                                        vecAbsoluteTolerance,
                                                        vecRelativeTolerance)
            satPropagator = NumericalPropagator(thisintegrator)
            satPropagator.setInitialState(SpacecraftState(initialCartesianOrbit, mass))
            satPropagator.setAttitudeProvider(nadirPointing)

        # determine the level of the propagator
        try:
            if (level == "Low") or (level == "low") or (level == "L") or (level == "l"):
                propagatorCase = 1
            elif (level == "Medium") or (level == "medium") or (level == "M") or (level == "m"):
                propagatorCase = 2
            elif (level == "High") or (level == "high") or (level == "H") or (level == "h"):
                propagatorCase = 3
            else:
                propagatorCase = 0
        except:
            propagatorCase = 0

        if (propagatorCase == 0):
            gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(1, 1)
            gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecefFrame, gravityProvider)

            satPropagator.addForceModel(gravityAttractionModel)


        elif (propagatorCase == 1):
            gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(2, 2)
            gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecefFrame, gravityProvider)

            satPropagator.addForceModel(gravityAttractionModel)

        elif (propagatorCase == 2):
            gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(21, 21)
            gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecefFrame, gravityProvider)
            # Atmospheric drag
            msafe = MarshallSolarActivityFutureEstimation(
                MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES,
                MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE)
            atmosphere = DTM2000(msafe, sun, wgs84Ellipsoid)
            isotropicDrag = IsotropicDrag(0.02, 2.2)
            dragForce = DragForce(atmosphere, isotropicDrag)
            # Solar radiation pressure
            isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(0.02, 1.0)
            solarRadiationPressure = SolarRadiationPressure(sun, wgs84Ellipsoid.getEquatorialRadius(),
                                                            isotropicRadiationSingleCoeff)
            # Third body attraction model
            moon_3dbodyattraction = ThirdBodyAttraction(moon)
            sun_3dbodyattraction = ThirdBodyAttraction(sun)

            satPropagator.addForceModel(gravityAttractionModel)
            satPropagator.addForceModel(dragForce)
            satPropagator.addForceModel(solarRadiationPressure)
            satPropagator.addForceModel(moon_3dbodyattraction)
            satPropagator.addForceModel(sun_3dbodyattraction)


        elif (propagatorCase == 3):

            # Earth gravity field with degree 64 and order 64
            gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(64, 64, self.startTime)
            gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecefFrame, gravityProvider)

            # Third body attraction model
            moon_3dbodyattraction = ThirdBodyAttraction(moon)
            sun_3dbodyattraction = ThirdBodyAttraction(sun)

            # Solar radiation pressure
            isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(0.02, 1.0)
            solarRadiationPressure = SolarRadiationPressure(sun, wgs84Ellipsoid, isotropicRadiationSingleCoeff)

            # Relativity
            relativity = Relativity(Constants.EIGEN5C_EARTH_MU)

            oceanicTides = OceanTides(FramesFactory.getITRF(IERSConventions.IERS_2010, True),
                                      Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_MU,
                                      5, 5, IERSConventions.IERS_2010,
                                      TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True))
            solidTidess = SolidTides(FramesFactory.getITRF(IERSConventions.IERS_2010, True),
                                     Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_MU,
                                     gravityProvider.getTideSystem(),
                                     IERSConventions.IERS_2010,
                                     TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True),
                                     [CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()])

            # Atmospheric drag
            # from org.orekit.models.earth.atmosphere import NRLMSISE00
            # atmosphere = NRLMSISE00(msafe, sun, wgs84Ellipsoid)
            msafe = MarshallSolarActivityFutureEstimation(
                MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES,
                MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE)
            atmosphere = DTM2000(msafe, sun, wgs84Ellipsoid)
            isotropicDrag = IsotropicDrag(0.02, 2.2)
            dragForce = DragForce(atmosphere, isotropicDrag)

            satPropagator.addForceModel(gravityAttractionModel)
            satPropagator.addForceModel(moon_3dbodyattraction)
            satPropagator.addForceModel(sun_3dbodyattraction)
            satPropagator.addForceModel(solarRadiationPressure)
            satPropagator.addForceModel(relativity)
            satPropagator.addForceModel(dragForce)
            satPropagator.addForceModel(oceanicTides)
            satPropagator.addForceModel(solidTidess)

        return satPropagator
    

    def rearangeTexture(self, texture):
        N_lat = int(texture.shape[0])
        N_lon = int(texture.shape[1])
        newTexture = np.zeros([N_lat, N_lon], dtype=np.uint8)
        newTexture[:1024] = texture[1024:]
        newTexture[1024:] = texture[:1024]
        return newTexture


    def sphereGenerator(self, size, texture): 
        N_lat = int(texture.shape[0])
        N_lon = int(texture.shape[1])
        theta = np.linspace(0,2*np.pi,N_lat)
        phi = np.linspace(0,np.pi,N_lon)
        
        # Set up coordinates for points on the sphere
        x0 = size * np.outer(np.cos(theta),np.sin(phi))
        y0 = size * np.outer(np.sin(theta),np.sin(phi))
        z0 = size * np.outer(np.ones(N_lat),np.cos(phi))
        
        # Set up trace
        return x0,y0,z0
    
    def getBlueMarbleFigure(self, path):
        earthMapTextureT = np.asarray(Image.open(path.format('earth'))).T
        earthMapTexture = self.rearangeTexture(earthMapTextureT)
        

        earthSphereX,earthSphereY,earthSphereZ =\
              self.sphereGenerator(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,earthMapTexture)
        earhtMarbleSurf = go.Surface(x=earthSphereX, y=earthSphereY, z=earthSphereZ,
                  surfacecolor=earthMapTexture,
                  colorscale= mainScenario.earthMapColorScale) 

        earhtMarbleFigure = go.Figure(data=[earhtMarbleSurf])

        return earhtMarbleFigure
    

    def ExBuildWalker( num_plane, num_sat, F, refSat):
        # allsat = [[0 for i in range(num_sat)] for i in range(num_plane)]
        orbits = []
        raan0 = refSat.getRightAscensionOfAscendingNode() * ( 180.0 / np.pi )
        ta0 = refSat.getTrueAnomaly() * ( 180.0 / np.pi )
        for i in range(num_plane):
            for j in range(num_sat):
                raan = raan0 + i * 360.0 / num_plane
                ta = ta0 + j * 360.0 / num_sat + i * 360 * F / (num_sat * num_plane)
                ta = ta % 360.0
                if ta >= 180.0:
                    ta -= 360.0
                if raan > 360.0:
                    raan -= 360.0
                newOrbit = KeplerianOrbit(refSat.getA(), refSat.getE(), refSat.getI(), 
                                          refSat.getPerigeeArgument(), raan * (np.pi / 180.0), 
                                          ta * ( np.pi / 180 ), PositionAngleType.TRUE,
                                           refSat.getFrame(), refSat.getPVCoordinates().getDate(), 
                                           Constants.WGS84_EARTH_MU)

                orbits.append(newOrbit)
        return orbits

    def buildISLStar(self, num_plane, num_sat, F, activeSats):
        
        '''
        ** reference satellite must be assigned before calling this method!!!
        inputs: number of plane: int
                number of satellites per plane: int
                walker phasing parameter: int (between zero and np-1)
                reference satellite: orekit KeplerianOrbit
                list of satellites with ISL link: list[dtype=int] length(between 1 and np*ns)
        '''

        if (num_plane-np.floor(num_plane)!=0):
            num_plane = np.floor(num_plane)

        if (num_sat-np.floor(num_sat)!=0):
            num_sat = np.floor(num_sat)
        
        if (F-np.floor(F)!=0):
            F = np.floor(F)
        
        if F < 0 :
            F = 0

        if F > num_plane:
            F = num_plane
        

        for i in range(len(activeSats)):
            if (activeSats[i]-np.floor(activeSats[i])!=0 ):
                activeSats[i]=np.floor(activeSats[i])
            if  activeSats[i] < 0:
                 activeSats[i] = 0

            if  activeSats[i] > num_plane*num_sat:
                 activeSats[i] = num_plane*num_sat   

        activeSats = np.unique((activeSats))

        nominalLinksList = np.zeros([len(activeSats),5],dtype=int)
        allContributedSatellitesList = []
        selectedOrbitsIndexes = []
        selectedOrbits = []
        orbits = []

        refSat = self.referenceOrbit
        raan0 = refSat.getRightAscensionOfAscendingNode() * (180.0 / np.pi)
        ta0 = refSat.getTrueAnomaly() * (180.0 / np.pi)
        for i in range(num_plane):
            for j in range(num_sat):
                raan = raan0 + i * 180.0 / num_plane
                ta = ta0 + j * 360.0 / num_sat + i * 360 * F / (num_sat * num_plane)
                ta = ta % 360.0

                if raan > 180.0:
                    raan -= 180.0
                newOrbit = KeplerianOrbit(refSat.getA(), refSat.getE(), refSat.getI(),
                                          refSat.getPerigeeArgument(), raan * (np.pi / 180.0),
                                          ta * (np.pi / 180), PositionAngleType.TRUE,
                                          refSat.getFrame(), refSat.getPVCoordinates().getDate(),
                                          Constants.WGS84_EARTH_MU)
                orbits.append(newOrbit)
        print(" ")
        print("walker-star constellation created ({}, {}, {}).".format(num_sat*num_plane, num_plane, F))

        for i in range(len(activeSats)):
            currentSatellite = activeSats[i]
            nominalLinksList[i , 0] = currentSatellite
            # finding neigbhor satellites 
            currentPlane = int(np.floor((currentSatellite - 1)/(num_sat) ) + 1)

            if (currentPlane == 1):
                nominalLinksList[i, 1] = -1
                nominalLinksList[i, 3] = currentSatellite + num_sat
            elif (currentPlane == num_plane):
                nominalLinksList[i, 1] = currentSatellite - num_sat
                nominalLinksList[i, 3] = -1
            else: 
                nominalLinksList[i, 1] = currentSatellite - num_sat
                nominalLinksList[i, 3] = currentSatellite + num_sat
            
            if  (currentSatellite == (currentPlane - 1)*num_sat + 1):
                nominalLinksList[i, 2] = currentSatellite + 1
                nominalLinksList[i, 4] = currentSatellite + num_sat - 1
            elif (currentSatellite == currentPlane*num_sat):
                nominalLinksList[i, 2] = currentSatellite - num_sat + 1
                nominalLinksList[i, 4] = currentSatellite - 1
            else:
                nominalLinksList[i, 2] = currentSatellite + 1
                nominalLinksList[i, 4] = currentSatellite - 1

            print(" ")
            print("Satellite {} is equipped with ISL links: {}".format(currentSatellite , nominalLinksList[i , 1:]))
            print("         {}         ".format(nominalLinksList[i, 2]))
            print("          |          ")
            print("          |          ")
            print("{}-------{}---------{}".format(nominalLinksList[i, 1], nominalLinksList[i, 0], nominalLinksList[i, 3]))
            print("          |          ")
            print("          |          ")
            print("         {}          ".format(nominalLinksList[i, 4]))
            
        for i in range(len(activeSats)):
            allContributedSatellitesList.append(activeSats[i])
            for j in range(4):
                if nominalLinksList[i , j + 1] == -1:
                    pass
                else:
                    allContributedSatellitesList.append(nominalLinksList[i, j + 1])
        uniqueIndexes = np.unique(np.array(allContributedSatellitesList))
        for i in range(len(uniqueIndexes)):
            selectedOrbits.append(orbits[uniqueIndexes[i] - 1])
            selectedOrbitsIndexes.append(uniqueIndexes[i])

        return orbits, selectedOrbits, selectedOrbitsIndexes, nominalLinksList

class MyObserver(PythonBatchLSObserver):
    def evaluationPerformed(self, itCounts,
                            evCounts, orbits, orbParams, propParams,
                            measParams, provider, lspEval):
        print(itCounts)

# building scenario
mainScenario = ssEnvironmentBuilder()
duration = mainScenario.duration
mainScenario.addReferenceOrbit(700000.0, 0.0, 89.46, 0.0, 0.0, 0.0)
#mainScenario.feedGroundStations(["Tabriz", "Boushehr", "Bandar Abbas", "Mashhad"])

underODSatellites = [100]
starConstellation, activeSatellites, selectedOrbitsIndexes, linkList = \
    mainScenario.buildISLStar( 15, 28, 1, underODSatellites)

if mainScenario.showActiveSatellites:
    mainScenario.displayActiveSatellitesConfiguration(activeSatellites, selectedOrbitsIndexes, linkList)

underODOrbits = []
linkedOrbits = []
allOrbits = []

for i in range(len(underODSatellites)):
    underODOrbits.append(starConstellation[underODSatellites[i] - 1])
    allOrbits.append(starConstellation[underODSatellites[i] - 1])
    linkedTemp = []
    for j in range(len(linkList[i])):
        if (linkList[i, j] != -1) and (linkList[i, j] != underODSatellites[i]):
            linkedTemp.append(starConstellation[int(linkList[i, j] - 1)])
            allOrbits.append(starConstellation[int(linkList[i, j] - 1)])
        else:
            pass
    linkedOrbits.append(linkedTemp)

underODPropagators = []
linkedPropagators = []
for i in range(len(underODOrbits)):
    underODPropagators.append(mainScenario.getPropagator("Get", underODOrbits[i], 100.0, "H"))
    linkedTemp = []
    lt = linkedOrbits[i]
    for j in range(len(lt)):
        linkedTemp.append(mainScenario.getPropagator("Get", lt[j], 100.0, "H"))
    linkedPropagators.append(linkedTemp) 

propBuilder = []
propBuilder.append(mainScenario.getPropagator("Builder", underODOrbits[0], 100.0, "H"))

thisLinkList = linkList[0]
linked = linkedOrbits[0]
for i in range(len(linked)):
    propBuilder.append(mainScenario.getPropagator("Builder", linked[i], 100.0, "H"))

thisEstimator = mainScenario.getPodEstimator(propBuilder)
ODStartTime = underODOrbits[0].getDate()

linked = []
linked = linkedPropagators[0]

if mainScenario.isISL:
    mainScenario.feedIslMeasurements(thisEstimator, underODPropagators[0], linked, ODStartTime, duration, 'Range', thisLinkList)
    mainScenario.feedIslMeasurements(thisEstimator, underODPropagators[0], linked, ODStartTime, duration, 'Phase', thisLinkList)

mainScenario.feedPodMeasurements(thisEstimator,allOrbits, ODStartTime, duration, "GPS", thisLinkList)

estimatorOutput = thisEstimator.estimate()
estimatedOrbits = []
estimatedPropagators = []
realPropagatores = []
for indx in range(len(propBuilder)):
    thisOrbit = []
    thisOrbit = estimatorOutput[indx].getInitialState().getOrbit()
    estimatedOrbits.append(estimatorOutput[indx].getInitialState().getOrbit())
    midEstProp = mainScenario.getPropagator("Get", thisOrbit, 100.0, "H")
    midRealProp = mainScenario.getPropagator("Get", allOrbits[indx], 100.0, "H")
    estEphemerisGenerator = midEstProp.getEphemerisGenerator()
    realEphemerisGenerator = midRealProp.getEphemerisGenerator()
    midEstProp.propagate(mainScenario.startTime, mainScenario.endTime)
    midRealProp.propagate(mainScenario.startTime, mainScenario.endTime)
    estimatedPropagators.append(estEphemerisGenerator.getGeneratedEphemeris())
    realPropagatores.append(realEphemerisGenerator.getGeneratedEphemeris())
    ricposdiff = mainScenario.getRICPositionDifference(realPropagatores[indx], estimatedPropagators[indx],
                                                        mainScenario.startTime, mainScenario.endTime, mainScenario.plotStepTime)
    if mainScenario.isISL:
        additionalFigureInfo = 'ISL is ON, ISL sigma = ' + str(mainScenario.islSigma) + \
        ', ISL BaseWeight = ' + str(mainScenario.islBaseWeight) + ', GPS sigma = ' + str(mainScenario.gpsSigma) + \
        ', GPS BaseWeight = ' + str(mainScenario.gpsBaseWeight) + ', GPS noise sigma = ' + str(mainScenario.gpsPositionStd) + \
            ',  for Sat:#' +  str(thisLinkList[indx])
    else:
        additionalFigureInfo = 'ISL is OFF, GPS sigma = ' + str(mainScenario.gpsSigma) + \
        ', GPS BaseWeight = ' + str(mainScenario.gpsBaseWeight) + ' GPS noise sigma = ' + str(mainScenario.gpsPositionStd) + \
            ',  for Sat:#' +  str(thisLinkList[indx])
    resifig = mainScenario.getResidualFigure(ricposdiff, additionalFigureInfo)
    pio.show(resifig)