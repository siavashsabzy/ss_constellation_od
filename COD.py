import orekit

orekit.initVM()
from orekit import JArray, JArray_double
from org.orekit.models.earth import ReferenceEllipsoid
from org.orekit.attitudes import NadirPointing
from org.orekit.propagation.events import ElevationDetector
from org.orekit.propagation.events.handlers import ContinueOnEvent
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.gravity import ThirdBodyAttraction
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
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
from org.orekit.estimation.measurements import GroundStation, Range, RangeRate, ObservableSatellite, AngularAzEl, \
    Position
from org.orekit.estimation.measurements.gnss import OneWayGNSSRange
from org.orekit.estimation.measurements.generation import Generator, RangeBuilder, RangeRateBuilder, AngularAzElBuilder, \
    EventBasedScheduler, SignSemantic
# from org.orekit.estimation.measurements.generation import RangeBuilder, RangeRateBuilder
from org.orekit.frames import FramesFactory, TopocentricFrame, ITRFVersion, LOFType
from org.orekit.time import AbsoluteDate, TimeScalesFactory, FixedStepSelector
from org.orekit.utils import Constants, IERSConventions, PVCoordinates, TimeStampedPVCoordinates
from org.orekit.bodies import GeodeticPoint, CelestialBodyFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngle, CartesianOrbit
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
from org.hipparchus.linear import MatrixUtils

# Atmospheric Models
from org.orekit.models.earth import EarthITU453AtmosphereRefraction
# Tropospheric Models
from org.orekit.models.earth.troposphere import SaastamoinenModel

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
from org.orekit.propagation.analytical.gnss import GNSSPropagator, GNSSPropagatorBuilder
from org.orekit.propagation.analytical.gnss.data import CommonGnssData, GNSSOrbitalElements, AbstractNavigationMessage, GPSNavigationMessage
from org.orekit.propagation.conversion import PropagatorConverter, FiniteDifferencePropagatorConverter, AbstractPropagatorConverter


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


class GroundStation_:
    def __init__(self, name, latitude, longitude, altitude, elevationAngle, temperature, humidity, pressure,
                 identification):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.elevationAngle = (elevationAngle) * (np.pi / 180.0)
        self.temperature = temperature
        self.humidity = humidity
        self.pressure = pressure
        self.identification = identification

    def groundStation(self, ellipsoid):
        geodeticPoint = GeodeticPoint(self.latitude, self.longitude, self.altitude)
        topocentricFrame = TopocentricFrame(ellipsoid, geodeticPoint, self.name)
        return GroundStation(topocentricFrame)

    def atmosphericRefraction(self):
        return EarthITU453AtmosphereRefraction(self.altitude)

    def troposphericModel(self):
        return SaastamoinenModel(self.temperature, self.pressure, self.humidity)

# def ionosphericModel(self):


class COD_():

    def __init__(self):
        super().__init__()
        self.RAD2DEG = (180 / np.pi)
        self.DEG2RAD = (np.pi / 180)
        self.rangeMaxError = 0.001
        self.rangeSigma = 0.01
        self.rangeWeight = 0.01
        self.groundStationsList = []
        self.numberOfGroundStation = 0
        self.spacecraftsList = []
        self.isTwoWay = True
        self.ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
        self.eciFrame = FramesFactory.getGCRF()
        self.wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(self.ecefFrame)

    def inter(self):
        pass

    def addGroundStation(self, latitude, longitude, altitude, name):
        stationPoint = GeodeticPoint(latitude * self.DEG2RAD, longitude * self.DEG2RAD, altitude)
        stationFrame = TopocentricFrame(self.wgs84Ellipsoid, stationPoint, name)
        groundStation = GroundStation(stationFrame)
        self.groundStationsList.append(groundStation)
        self.numberOfGroundStation += 1

    def addSpacecraft(self, spacecraft):
        self.spacecraftsList.append(spacecraft)

    def generate_measurements(propagator, station, station_name, altitude, elevAngle, meas_type, t0, duration,
                              sigma=0.00400,  # RangeRate sigma
                              base_weight=1.0,
                              two_way=True,
                              withRefraction=True,
                              step=1.0,
                              seed=0.0):

        measurementslist = []
        tf = t0.shiftedBy(3600.0 * duration)

        small = 1e-10
        random_generator = Well19937a(int(seed))
        gaussian_generator = GaussianRandomGenerator(random_generator)
        covariance = MatrixUtils.createRealDiagonalMatrix([float(sigma * sigma), float(sigma * sigma)])
        noise_source = CorrelatedRandomVectorGenerator(covariance, float(small), gaussian_generator)

        generator = Generator()
        satellite = generator.addPropagator(propagator)

        if meas_type == 'RANGE':
            builder = RangeBuilder(noise_source, station, two_way, sigma, base_weight, satellite)
        elif meas_type == 'RANGERATE':
            builder = RangeRateBuilder(noise_source, station, two_way, sigma, base_weight, satellite)
        elif meas_type == 'AZEL':
            sigmaAzEl = JArray_double([sigma, sigma])
            weightAzEl = JArray_double([base_weight, base_weight])
            builder = AngularAzElBuilder(noise_source, station, sigmaAzEl, weightAzEl, satellite)

        fixed_step_selector = FixedStepSelector(step, TimeScalesFactory.getUTC())
        if withRefraction:
            elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(
                elevAngle).withRefraction(EarthITU453AtmosphereRefraction(altitude)).withHandler(ContinueOnEvent())
        else:
            elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(elevAngle).withHandler(
                ContinueOnEvent())
        scheduler = EventBasedScheduler(builder, fixed_step_selector, propagator, elevation_detector,
                                        SignSemantic.FEASIBLE_MEASUREMENT_WHEN_POSITIVE)
        generator.addScheduler(scheduler)

        measurements = generator.generate(t0, tf)
        for measObject in measurements:
            if meas_type == 'RANGE':
                meas = Range.cast_(measObject)
            elif meas_type == 'RANGERATE':
                meas = RangeRate.cast_(measObject)
            elif meas_type == 'AZEL':
                meas = AngularAzEl.cast_(measObject)
            else:
                raise ValueError("unrecognized measurement type")
            measurementslist.append(meas)
            if two_way:
                way_str = "TWOWAY"
            else:
                way_str = "ONEWAY"
            print("{}   {} {}       {}        {}".format(meas.getDate().toString(), way_str, meas_type, station_name,
                                                         meas.getObservedValue()[0]))
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

    def ExGetPropagator(propagatorType,
                        initialOrbit, mass,
                        eciFrame, ecefFrame,
                        level, minStep, maxStep,
                        vecAbsoluteTolerance, vecRelativeTolerance):

        moon = CelestialBodyFactory.getMoon()
        sun = CelestialBodyFactory.getSun()
        wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecefFrame)
        nadirPointing = NadirPointing(eciFrame, wgs84Ellipsoid)

        initialCartesianOrbit = CartesianOrbit(SpacecraftState(initialOrbit, mass).getPVCoordinates(eciFrame),
                                               eciFrame, wgs84Ellipsoid.getGM())

        if (propagatorType == 'Builder'):
            integratorBuilder = DormandPrince853IntegratorBuilder(minStep, maxStep, vecAbsoluteTolerance)

            satPropagator = NumericalPropagatorBuilder(initialCartesianOrbit, integratorBuilder, PositionAngle.TRUE,
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
            gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(64, 64)
            gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecefFrame, gravityProvider)

            # Third body attraction model
            moon_3dbodyattraction = ThirdBodyAttraction(moon)
            sun_3dbodyattraction = ThirdBodyAttraction(sun)

            # Solar radiation pressure
            isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(0.02, 1.0)
            solarRadiationPressure = SolarRadiationPressure(sun, wgs84Ellipsoid.getEquatorialRadius(),
                                                            isotropicRadiationSingleCoeff)

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

    def ExBuildWalker( num_plane, num_sat, F, refSat, activeSats):
        '''
        inputs: number of plane: int
                number of satellites per plane: int
                walker phasing parameter: int (between zero and np-1)
                reference satellite: orekit KeplerianOrbit
                list of satellites with ISL link: list[dtype=int] length(between 1 and np*ns)
        '''
        # allsat = [[0 for i in range(num_sat)] for i in range(num_plane)]
        nominalLinksList = np.zeros([len(activeSats),4],dtype=int)
        allContributedSatellitesList = []
        selectedOrbitsIndexes = []
        selectedOrbits = []
        orbits = []
        raan0 = refSat.getRightAscensionOfAscendingNode() * (180.0 / np.pi)
        ta0 = refSat.getTrueAnomaly() * (180.0 / np.pi)
        for i in range(num_plane):
            for j in range(num_sat):
                raan = raan0 + i * 180.0 / num_plane
                ta = ta0 + j * 360.0 / num_sat + i * 360 * F / (num_sat * num_plane)
                ta = ta % 360.0
                #if ta >= 180.0:
                #    ta -= 360.0
                if raan > 180.0:
                    raan -= 180.0
                newOrbit = KeplerianOrbit(refSat.getA(), refSat.getE(), refSat.getI(),
                                          refSat.getPerigeeArgument(), raan * (np.pi / 180.0),
                                          ta * (np.pi / 180), PositionAngle.TRUE,
                                          refSat.getFrame(), refSat.getPVCoordinates().getDate(),
                                          Constants.WGS84_EARTH_MU)
                
                orbits.append(newOrbit)
        
        for i in range(len(activeSats)):
            currentSatellite = activeSats[i]
            
            # finding neigbhor satellites 
            currentPlane = np.floor(currentSatellite/num_sat) + 1

            if (currentPlane == 1):
                nominalLinksList[i, 3] = -1
                neighbourPlanes = [ currentPlane + 1, -1] 
                nominalLinksList[i, 2] = currentSatellite + num_sat
            elif (currentPlane == num_plane):
                nominalLinksList[i, 3] = -1
                neighbourPlanes = [ currentPlane - 1 , -1]
                nominalLinksList[i, 2] = currentSatellite - num_sat
            else: 
                neighbourPlanes = [ currentPlane - 1 , currentPlane + 1]
                nominalLinksList[i, 2] = currentSatellite - num_sat
                nominalLinksList[i, 3] = currentSatellite + num_sat



            if  (currentSatellite == (currentPlane - 1)*num_sat + 1):
                nominalLinksList[i, 0] = currentSatellite + 1
                nominalLinksList[i, 1] = currentSatellite + num_sat - 1
            elif (currentSatellite == currentPlane*num_sat):
                nominalLinksList[i, 0] = currentSatellite - 1
                nominalLinksList[i, 1] = currentSatellite - num_sat + 1
            else:
                nominalLinksList[i, 0] = currentSatellite - 1
                nominalLinksList[i, 1] = currentSatellite + 1
            
        for i in range(len(activeSats)):
            allContributedSatellitesList.append(activeSats[i])
            for j in range(4):
                if nominalLinksList[i, j] == -1:
                    pass
                else:
                    allContributedSatellitesList.append(nominalLinksList[i, j])
        uniqueIndexes = np.unique(np.array(allContributedSatellitesList))
        for i in range(len(uniqueIndexes)):
            selectedOrbits.append(orbits[uniqueIndexes[i] - 1])
            selectedOrbitsIndexes.append(uniqueIndexes[i])
        
        return orbits, selectedOrbits,selectedOrbitsIndexes, nominalLinksList


class MyObserver(PythonBatchLSObserver):
    def evaluationPerformed(self, itCounts,
                            evCounts, orbits, orbParams, propParams,
                            measParams, provider, lspEval):
        print(itCounts)


"""
# ###################################################################################
# ###################################################################################
##########  Test Area ##############**********************************


utc = TimeScalesFactory.getUTC()
#startTime = AbsoluteDate(2022, 10, 3, 8, 30, 0.0, utc)
eciFrame = FramesFactory.getGCRF()
ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)



#refsat_test = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 500000.0,
#                0.0 , 78.0*np.pi/180, 
#                0.0*np.pi/180,  180.0*np.pi/180 ,  25.0*np.pi/180,  PositionAngle.TRUE,
#                eciFrame,  startTime, Constants.WGS84_EARTH_MU)

#cons = COD_().ExBuildWalker(25, 13, 1, refsat_test)

with open('S1A_OPER_AUX_RESORB_OPOD_20220103T002736_V20220102T202842_20220102T234612.EOF', 'r') as ephfile:
    sentineldata = ephfile.read()

Bs_sentineldata = BeautifulSoup(sentineldata, "xml")
b_unique = Bs_sentineldata.find_all('OSV')
position = pd.DataFrame(columns=['X', 'Y', 'Z','Vx','Vy','Vz'])

for i in range(len(b_unique)):
    s = b_unique[i].find_all("UTC")[0].contents[0]
    x = b_unique[i].find_all("X")[0].contents[0]
    y = b_unique[i].find_all("Y")[0].contents[0]
    z = b_unique[i].find_all("Z")[0].contents[0]
    vx = b_unique[i].find_all("VX")[0].contents[0]
    vy = b_unique[i].find_all("VY")[0].contents[0]
    vz = b_unique[i].find_all("VZ")[0].contents[0]
    a = datetime.strptime(s, "UTC=%Y-%m-%dT%H:%M:%S.%f")
    position.loc[a] = [float(x), float(y), float(z), float(vx), float(vy), float(vz)]

d = position.iloc[0].name



thisPV = PVCoordinates(Vector3D(5665849.027545, -4240370.657895,-40753.628053),
                                   Vector3D(-921.855016,-1288.326693,7430.094083))
thisPVT = TimeStampedPVCoordinates(AbsoluteDate(2022, 1, 2, 20, 28, 42.382456, utc),
                                   Vector3D(5665849.027545, -4240370.657895,-40753.628053),
                                   Vector3D(-921.855016,-1288.326693,7430.094083))
ecef2eci = ecefFrame.getTransformTo(eciFrame,AbsoluteDate(2022, 1, 2, 20, 28, 42.382456, utc))
inertialPV = ecef2eci.transformPVCoordinates(thisPV)
inertialPVT = TimeStampedPVCoordinates(AbsoluteDate(2022, 1, 2, 20, 28, 42.382456, utc),
                                        inertialPV.getPosition(),
                                        inertialPV.getVelocity())

thisOrbit = CartesianOrbit(inertialPV, eciFrame,
                            AbsoluteDate(2022, 1, 2, 20, 28, 42.382456, utc),
                           Constants.WGS84_EARTH_MU )
sc = SpacecraftState(thisOrbit)

prop = COD_().ExGetPropagator(sc, "H", 0.1,1.0,0.01,0.01)

secpv = prop.propagate(AbsoluteDate(2022, 1, 2, 23, 46, 12.382456, utc)).getPVCoordinates()

eci2ecef = eciFrame.getTransformTo(ecefFrame,AbsoluteDate(2022, 1, 2, 23, 46, 12.382456, utc))
resu = eci2ecef.transformPVCoordinates(secpv)
print(resu.getPosition())
print(resu.getVelocity())
a =0

#print(cons[1].getRightAscensionOfAscendingNode())



###########################################################################
############################################################################
############################################################################
"""

# initial figure configuration
earthMapColorScale =[[0.0, 'rgb(30, 59, 117)'],
                 [0.1, 'rgb(46, 68, 21)'],
                 [0.2, 'rgb(74, 96, 28)'],
                 [0.3, 'rgb(115,141,90)'],
                 [0.4, 'rgb(122, 126, 75)'],
                 [0.6, 'rgb(122, 126, 75)'],
                 [0.7, 'rgb(141,115,96)'],
                 [0.8, 'rgb(223, 197, 170)'],
                 [0.9, 'rgb(237,214,183)'],
                 [1.0, 'rgb(255, 255, 255)']]

def sphereGenerator(size, texture): 
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

def rearangeTexture(texture):
    N_lat = int(texture.shape[0])
    N_lon = int(texture.shape[1])
    newTexture = np.zeros([N_lat, N_lon], dtype=np.uint8)
    newTexture[:1024] = texture[1024:]
    newTexture[1024:] = texture[:1024]
    return newTexture
    


earthMapTextureT = np.asarray(Image.open('earth.jpeg'.format('earth'))).T
earthMapTexture = rearangeTexture(earthMapTextureT)


earthSphereX,earthSphereY,earthSphereZ = sphereGenerator(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,earthMapTexture)
earhtMarbleSurf = go.Surface(x=earthSphereX, y=earthSphereY, z=earthSphereZ,
                  surfacecolor=earthMapTexture,
                  colorscale=earthMapColorScale) 

earhtMarbleFigure = go.Figure(data=[earhtMarbleSurf])







elevationMask = 10.0  # in degrees

# Scenario Start Time
utc = TimeScalesFactory.getUTC()
startTime = AbsoluteDate(2022, 10, 3, 8, 30, 0.0, utc)
endTime = AbsoluteDate(2022, 10, 3, 8, 50, 0.0, utc)
# referenceDate = AbsoluteDate.J2000_EPOCH  # reference date
# mjdUtcEpoch = AbsoluteDate(1858, 11, 17, 0, 0, 0.0, utc)


ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
eciFrame = FramesFactory.getGCRF()
moon = CelestialBodyFactory.getMoon()
sun = CelestialBodyFactory.getSun()
wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecefFrame)



#stationPoint = GeodeticPoint(35.0 * np.pi / 180.0, 51.0 * np.pi / 180.0, 100.0)
#ecefpos = TopocentricFrame(wgs84Ellipsoid, stationPoint, 'name').getPVCoordinates(startTime, ecefFrame).getPosition()
#earhtMarbleFigure.add_trace(go.Scatter3d(x = [ecefpos.getX()], 
#                                         y = [ecefpos.getY()], 
#                                         z = [ecefpos.getZ()]))
#pio.show(earhtMarbleFigure)

initialRealOrbitSat1 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 700000.0,  # Semimajor Axis (m)
                                      0.0,  # Eccentricity
                                      89.0 * np.pi / 180,  # Inclination (radians)
                                      0.0 * np.pi / 180,  # Perigee argument (radians)
                                      0.0 * np.pi / 180,  # Right ascension of ascending node (radians)
                                      0.0 * np.pi / 180,  # Anomaly (rad/s)
                                      PositionAngle.TRUE,  # Sets which type of anomaly we use
                                      eciFrame,
                                      # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                                      startTime,  # Sets the date of the orbital parameters
                                      Constants.WGS84_EARTH_MU)  # Sets the central attraction coefficient (m³/s²)
initialRealOrbitSat2 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 700000.0,  # Semimajor Axis (m)
                                      0.0,  # Eccentricity
                                      78.0 * np.pi / 180,  # Inclination (radians)
                                      0.0 * np.pi / 180,  # Perigee argument (radians)
                                      207.0 * np.pi / 180,  # Right ascension of ascending node (radians)
                                      25.0 * np.pi / 180,  # Anomaly (rad/s)
                                      PositionAngle.TRUE,  # Sets which type of anomaly we use
                                      eciFrame,
                                      # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                                      startTime,  # Sets the date of the orbital parameters
                                      Constants.WGS84_EARTH_MU)  # Sets the central attraction coefficient (m³/s²)
initialEstimatedOrbitSat1 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 501234.0,
                                           0.0, 78.1 * np.pi / 180, 0.0 * np.pi / 180, 179.5 * np.pi / 180,
                                           26.00 * np.pi / 180,
                                           PositionAngle.TRUE, eciFrame, startTime, Constants.WGS84_EARTH_MU)
initialEstimatedOrbitSat2 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 500221.0,
                                           0.0, 78.0 * np.pi / 180, 0.0 * np.pi / 180, 207.0 * np.pi / 180,
                                           25.0 * np.pi / 180,
                                           PositionAngle.TRUE, eciFrame, startTime, Constants.WGS84_EARTH_MU)


starConstellation, activeSatellites, selectedOrbitsIndexes, linkList = COD_.ExBuildWalker( 15, 28, 1, initialRealOrbitSat1, [30])
epochECI2ECEF = eciFrame.getTransformTo(ecefFrame, startTime)
initialConfig = []
for j in range(len(activeSatellites)):
    eciPos = activeSatellites[j].getPVCoordinates().getPosition()
    truePositionECEF = epochECI2ECEF.getRotation().applyTo(eciPos)
    earhtMarbleFigure.add_trace(go.Scatter3d(x = [truePositionECEF.getX()], 
                                         y = [truePositionECEF.getY()], 
                                         z = [truePositionECEF.getZ()]))
    initialConfig.append(truePositionECEF)

earhtMarbleFigure.update_coloraxes(showscale=False)
earhtMarbleFigure.update_layout(showlegend=False)
pio.show(earhtMarbleFigure)

# Orbit propagator parameters
prop_min_step = 0.01  # s
prop_max_step = 2.0  # s
prop_position_error = 3.0  # m

sat1RealPropagator = COD_.ExGetPropagator("Get", initialRealOrbitSat1, 100.0,
                                          eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                          prop_position_error, prop_position_error)
sat2RealPropagator = COD_.ExGetPropagator("Get", initialRealOrbitSat2, 100.0,
                                          eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                          prop_position_error, prop_position_error)




sat1EstimatedPropagatorBuilder = COD_.ExGetPropagator("Builder", initialEstimatedOrbitSat1, 100.0,
                                                      eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                                      prop_position_error, prop_position_error)

# aa = FiniteDifferencePropagatorConverter(sat1EstimatedPropagatorBuilder, 1e-3, 1000)
# aa.convert()

sat2EstimatedPropagatorBuilder = COD_.ExGetPropagator("Builder", initialEstimatedOrbitSat2, 100.0,
                                                      eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                                      prop_position_error, prop_position_error)

estimator_convergence_thres = 1e-5
estimator_max_iterations = 25
estimator_max_evaluations = 35

matrixDecomposer = QRDecomposer(1e-11)
optimizer = GaussNewtonOptimizer(matrixDecomposer, False)
sat1Estimator = BatchLSEstimator(optimizer, sat1EstimatedPropagatorBuilder)
sat1Estimator.setParametersConvergenceThreshold(estimator_convergence_thres)
sat1Estimator.setMaxIterations(estimator_max_iterations)
sat1Estimator.setMaxEvaluations(estimator_max_evaluations)

##########################################################
# Station Filter

range_weight = 0.01  # Will be normalized later (i.e divided by the number of observations)
range_sigma = 0.01  # Estimated covariance of the range measurements, in meters

isTwoWay = True
currentDateTime = startTime

this_cons = COD_()
this_cons.addGroundStation(30.0, 48.5, 0.0, "Boushehr")
this_cons.addGroundStation(35.0, 51.0, 0.0, "Tehran")
#### this_cons.addGroundStation(30.0, 50.0, 0.0, "Boushehr")
#### this_cons.addGroundStation(33.0, 48.0, 1500.0, "Aleshtar")
#this_cons.addGroundStation(25.0, 59.0, 0.0, "Bandar Abbas")
#this_cons.addGroundStation(35.0, 58.0, 0.0, "Mashhad")
#this_cons.addGroundStation(38.0, 46.0, 0.0, "Tabriz")
#this_cons.addGroundStation(27.0, 63.0, 0.0, "Sistan")


for i in range(this_cons.numberOfGroundStation):
    meas = COD_.generate_measurements(sat1RealPropagator, this_cons.groundStationsList[i], str(i), 0.0,
                                      5.0 * np.pi / 180.0, "RANGERATE",
                                      currentDateTime, 1.0)
    for m in meas:
        sat1Estimator.addMeasurement(m)
"""
while currentDateTime.compareTo(endTime) <= 0:
    meas = this_cons.getRangeMeasurements(sat1RealPropagator, currentDateTime, elevationMask)
    currentDateTime = currentDateTime.shiftedBy(1.0)
    for m in meas:
        sat1Estimator.addMeasurement(m)
"""

sat1Estimator.setObserver(MyObserver())
print("Measurements feeded")
sat1EstimatedPropagator = sat1Estimator.estimate()[0]
sat1EstimatedInitialState = sat1EstimatedPropagator.getInitialState()
sat1EstimatedPropagator.resetInitialState(sat1EstimatedInitialState)
sat1EphemerisGenerator = sat1EstimatedPropagator.getEphemerisGenerator()

# Propagating from 1 day before data collection
# To 1 week after orbit determination (for CPF generation)
sat1EstimatedPropagator.propagate(startTime.shiftedBy(0.0), startTime.shiftedBy(2 * 3600.0 + 300))
sat1BoundedPropagator = sat1EphemerisGenerator.getGeneratedEphemeris()
initialState = sat1BoundedPropagator.getInitialState()

'''
brdctFirstGuess = AbstractNavigationMessage(Constants.WGS84_EARTH_MU, Constants.WGS84_EARTH_ANGULAR_VELOCITY, 10)
brdctFirstGuess.setSma(initialState.getA())
brdctFirstGuess.setSma(initialState.getE())
brdctFirstGuess.setI0(initialState.getI())
brdctFirstGuess.setTime(0.0)
brdctFirstGuess.setDate(initialState.getDate())
brdctFirstGuess.setPRN(1)
brdctFirstGuess.setM0(initialState.getKeplerianMeanMotion())

a = GNSSPropagatorBuilder(brdctFirstGuess)
'''


currentDateTime = startTime.shiftedBy(0.0)
endTime = startTime.shiftedBy(2 * 3600.0 + 300.0)

position_norm_resi = pd.DataFrame(columns=['PositionNormError'])
truePositionDataFrame = pd.DataFrame(columns=['x', 'y', 'z'])
truePositionECEFDataFrame = pd.DataFrame(columns=['x', 'y', 'z'])
estimatedPositionDataFrame = pd.DataFrame(columns=['x', 'y', 'z'])

while currentDateTime.compareTo(endTime) <= 0:
    thisECI2ECEF = eciFrame.getTransformTo(ecefFrame, currentDateTime)
    
    truePosition = sat1RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition()
    truePositionECEF = thisECI2ECEF.getRotation().applyTo(truePosition)
    estimatedPosition = sat1BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition()
    positionDifference = truePosition.subtract(estimatedPosition)
    errorInLVLH = LOFType.LVLH.rotationFromInertial(sat1RealPropagator.propagate(currentDateTime).getPVCoordinates()).applyTo(positionDifference)
    norm_resi = errorInLVLH.getNorm()

    truePositionDataFrame.loc[absolutedate_to_datetime(currentDateTime)] = \
        [truePosition.getX(), truePosition.getY(), truePosition.getZ()]
    truePositionECEFDataFrame.loc[absolutedate_to_datetime(currentDateTime)] = \
        [truePositionECEF.getX(), truePositionECEF.getY(), truePositionECEF.getZ()]
    estimatedPositionDataFrame.loc[absolutedate_to_datetime(currentDateTime)] = \
        [estimatedPosition.getX(), estimatedPosition.getY(), estimatedPosition.getZ()]
    position_norm_resi.loc[absolutedate_to_datetime(currentDateTime)] = norm_resi
    print(norm_resi)
    currentDateTime = currentDateTime.shiftedBy(10.0)

trace = go.Scattergl(
    x=position_norm_resi.index, y=position_norm_resi['PositionNormError'],
    mode='markers',
    name='PositionNormError'
)

truePositionTrace = go.Scatter3d( x = truePositionDataFrame['x'], 
                                  y = truePositionDataFrame['y'],
                                  z = truePositionDataFrame['z'])
estimatedPositionTrace = go.Scatter3d( x = estimatedPositionDataFrame['x'], 
                                        y = estimatedPositionDataFrame['y'],
                                        z = estimatedPositionDataFrame['z'])
truePositionECEFFigure = earhtMarbleFigure.add_trace(go.Scatter3d( x = truePositionECEFDataFrame['x'], 
                                                                y = truePositionECEFDataFrame['y'],
                                                                z = truePositionECEFDataFrame['z']))

data = [trace]

#tdLayout = go.Layout(title='3d trace',
#                      xaxis=dict(title='x [km]'), 
#                      yaxis=dict(title='y [km]'), 
#                      zaxis=dict(title='z [km]') )

layout = go.Layout(title='Position residuals', 
                   xaxis=dict(title='Datetime UTC'), 
                   yaxis=dict(title='Position residual (m)') )

fig = dict(data=data, layout=layout)
truePositionFigure = dict(data=truePositionTrace)
estimatedPositionFigure = dict(data=estimatedPositionTrace)
combinedFigure = dict(data=[truePositionTrace, estimatedPositionTrace])

# pio.write_image(fig, "file.png", height=1200, width=1600,scale=1)
pio.show(fig)
pio.show(truePositionECEFFigure)
pio.show(truePositionFigure)
pio.show(estimatedPositionFigure)
pio.show(combinedFigure)

sat2Estimator = BatchLSEstimator(optimizer, sat2EstimatedPropagatorBuilder)
sat2Estimator.setParametersConvergenceThreshold(estimator_convergence_thres)
sat2Estimator.setMaxIterations(estimator_max_iterations)
sat2Estimator.setMaxEvaluations(estimator_max_evaluations)
sat2Estimator.setObserver(MyObserver())
currentDateTime = startTime.shiftedBy(300.0)
endTime = startTime.shiftedBy(3600.0 + 300.0)
while currentDateTime.compareTo(endTime) <= 0:
    inter_range = sat2RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().distance(
        sat1RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition())
    oneway_meas = OneWayGNSSRange(sat1BoundedPropagator, 0.0,
                                  currentDateTime,
                                  inter_range,
                                  1.0,
                                  1.0,
                                  ObservableSatellite(0))
    currentDateTime = currentDateTime.shiftedBy(1.0)
    sat2Estimator.addMeasurement(oneway_meas)

sat2EstimatedPropagator = sat2Estimator.estimate()[0]
sat2EstimatedInitialState = sat2EstimatedPropagator.getInitialState()
sat2EstimatedPropagator.resetInitialState(sat2EstimatedInitialState)
sat2EphemerisGenerator = sat2EstimatedPropagator.getEphemerisGenerator()
sat2EstimatedPropagator.propagate(startTime.shiftedBy(300.0), startTime.shiftedBy(2 * 3600.0 + 300))
sat2BoundedPropagator = sat2EphemerisGenerator.getGeneratedEphemeris()
currentDateTime = startTime.shiftedBy(300.0)
endTime = startTime.shiftedBy(2 * 3600.0 + 300.0)
position_norm_resi = pd.DataFrame(columns=['PositionNormError'])
while currentDateTime.compareTo(endTime) <= 0:
    x_resi = (sat2RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getX() -
              sat2BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getX())
    y_resi = (sat2RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getY() -
              sat2BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getY())
    z_resi = (sat2RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getZ() -
              sat2BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getZ())
    norm_resi = (x_resi ** 2 + y_resi ** 2 + z_resi ** 2) ** 0.5
    position_norm_resi.loc[absolutedate_to_datetime(currentDateTime)] = norm_resi
    print(norm_resi)
    currentDateTime = currentDateTime.shiftedBy(10.0)

trace = go.Scattergl(
    x=position_norm_resi.index, y=position_norm_resi['PositionNormError'],
    mode='markers',
    name='PositionNormError'
)

data = [trace]

layout = go.Layout(
    title='Position residuals',
    xaxis=dict(
        title='Datetime UTC'
    ),
    yaxis=dict(
        title='Position residual (m)'
    )
)

fig = dict(data=data, layout=layout)

pio.show(fig)
# pio.write_image(fig, "test.svg", width=1.5*300, height=0.75*300, scale=1)
