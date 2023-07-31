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
from org.orekit.estimation.measurements import GroundStation, Range, RangeRate, ObservableSatellite, AngularAzEl, Position
from org.orekit.estimation.measurements.gnss import OneWayGNSSRange
from org.orekit.estimation.measurements.generation import Generator, RangeBuilder, RangeRateBuilder, AngularAzElBuilder, EventBasedScheduler, SignSemantic
#from org.orekit.estimation.measurements.generation import RangeBuilder, RangeRateBuilder
from org.orekit.frames import FramesFactory, TopocentricFrame, ITRFVersion
from org.orekit.time import AbsoluteDate, TimeScalesFactory, FixedStepSelector
from org.orekit.utils import Constants, IERSConventions, PVCoordinates, TimeStampedPVCoordinates
from org.orekit.bodies import GeodeticPoint, CelestialBodyFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngle, CartesianOrbit
from org.orekit.propagation import SpacecraftState
#from org.orekit.propagation.analytical import KeplerianPropagator
#from org.orekit.propagation.sampling import OrekitFixedStepHandler
from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation.numerical import NumericalPropagator
#from org.hipparchus.random import CorrelatedRandomVectorGenerator, NormalizedRandomGenerator
#from org.hipparchus.linear import RealMatrix
from org.hipparchus.geometry.euclidean.threed import Vector3D, SphericalCoordinates
from org.hipparchus.random import Well19937a, GaussianRandomGenerator, CorrelatedRandomVectorGenerator
from org.hipparchus.linear import MatrixUtils

# Atmospheric Models
from org.orekit.models.earth import EarthITU453AtmosphereRefraction
# Tropospheric Models
from org.orekit.models.earth.troposphere import SaastamoinenModel

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
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
webbrowser.register('firefox',None,webbrowser.BackgroundBrowser("C:\\Program Files\\Mozilla Firefox\\firefox.exe"))
pio.renderers.default = "firefox"


class GroundStation_:
	def __init__(self, name, latitude, longitude, altitude, elevationAngle, temperature, humidity, pressure, identification):
		self.name           = name
		self.latitude       = latitude
		self.longitude      = longitude
		self.altitude       = altitude
		self.elevationAngle = (elevationAngle) * (np.pi / 180.0)
		self.temperature    = temperature
		self.humidity       = humidity
		self.pressure       = pressure
		self.identification = identification

	def groundStation(self, ellipsoid):
		geodeticPoint    = GeodeticPoint(self.latitude, self.longitude, self.altitude)
		topocentricFrame = TopocentricFrame(ellipsoid, geodeticPoint, self.name)
		return GroundStation(topocentricFrame)

	def atmosphericRefraction(self):
		return EarthITU453AtmosphereRefraction(self.altitude)

	def troposphericModel(self):
		return SaastamoinenModel(self.temperature, self.pressure, self.humidity)

	#def ionosphericModel(self):




class COD_():

    def __init__(self):
        super().__init__()
        self.RAD2DEG = (180 / np.pi)
        self.DEG2RAD = (np.pi / 180)
        self.rangeMaxError = 0.001
        self.rangeSigma = 0.01
        self.rangeWeight = 0.01
        self.groundStationsList = []
        self.spacecraftsList = []
        self.isTwoWay = True
        self.ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
        self.eciFrame = FramesFactory.getGCRF()
        self.wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(self.ecefFrame)


    def inter(self):
        pass

    


    def addGroundStation(self, latitude, longitude, altitude, name):
        stationPoint  = GeodeticPoint(latitude*self.DEG2RAD, longitude*self.DEG2RAD, altitude)
        stationFrame  = TopocentricFrame(self.wgs84Ellipsoid, stationPoint, name)
        groundStation = GroundStation(stationFrame)
        self.groundStationsList.append(groundStation)

    def addSpacecraft(self, spacecraft):
        self.spacecraftsList.append(spacecraft)


    def generate_measurements(propagator, station, station_name, altitude, elevAngle, meas_type, t0, duration,
        sigma          = 0.040, # RangeRate sigma
        base_weight    = 1.0,
        two_way        = True,
        withRefraction = True,
        step           = 0.01,
        seed           = 0.0):

        measurementslist = []
        tf = t0.shiftedBy(3600.0*duration)
        
        small = 1e-10
        random_generator    = Well19937a(int(seed))
        gaussian_generator  = GaussianRandomGenerator(random_generator)
        covariance = MatrixUtils.createRealDiagonalMatrix([float(sigma * sigma),float(sigma * sigma)])
        noise_source = CorrelatedRandomVectorGenerator(covariance, float(small), gaussian_generator)

        generator  = Generator()
        satellite = generator.addPropagator(propagator)

        if meas_type == 'RANGE':
            builder = RangeBuilder(noise_source, station, two_way, sigma, base_weight, satellite)
        elif meas_type == 'RANGERATE':
            builder = RangeRateBuilder(noise_source, station, two_way, sigma, base_weight, satellite)
        elif meas_type == 'AZEL':
            sigmaAzEl  = JArray_double([sigma, sigma])
            weightAzEl = JArray_double([base_weight, base_weight])
            builder = AngularAzElBuilder(noise_source, station, sigmaAzEl, weightAzEl, satellite)

        fixed_step_selector = FixedStepSelector(step, TimeScalesFactory.getUTC())
        if withRefraction:
            elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(elevAngle).withRefraction(EarthITU453AtmosphereRefraction(altitude)).withHandler(ContinueOnEvent())
        else:
            elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(elevAngle).withHandler(ContinueOnEvent())
        scheduler = EventBasedScheduler(builder, fixed_step_selector, propagator, elevation_detector, SignSemantic.FEASIBLE_MEASUREMENT_WHEN_POSITIVE)
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
            print("{}   {} {}       {}        {}".format(meas.getDate().toString(), way_str, meas_type, station_name, meas.getObservedValue()[0]))
        return measurementslist


    def getRangeMeasurements(self, propagator, absoluteDateTime, elevationMask):
        currentInertialCoordinates = propagator.propagate(absoluteDateTime).getPVCoordinates()
        rangeMeasurementsList = []
        for groundStation in self.groundStationsList:
            elevation = groundStation.getBaseFrame().getElevation(
                currentInertialCoordinates.getPosition(), propagator.getFrame(), absoluteDateTime)
            if ((elevation * self.RAD2DEG) > elevationMask):
                realRange   = groundStation.getBaseFrame().getRange(currentInertialCoordinates.getPosition(), 
                                        propagator.getFrame(), absoluteDateTime)
                MeasuredRange = random.uniform(-self.rangeMaxError,self.rangeMaxError) + realRange

                Measurement = Range(groundStation, self.isTwoWay, absoluteDateTime, MeasuredRange,
                                self.rangeSigma, self.rangeWeight, ObservableSatellite(0) )
            
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
            
            satPropagator = NumericalPropagatorBuilder(initialCartesianOrbit, integratorBuilder, PositionAngle.TRUE, 1.0)
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
                                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,  Constants.WGS84_EARTH_MU,
                                        5, 5, IERSConventions.IERS_2010, 
                                    TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True))
            solidTidess = SolidTides(FramesFactory.getITRF(IERSConventions.IERS_2010, True),
                                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,  Constants.WGS84_EARTH_MU,
                                    gravityProvider.getTideSystem(),
                                    IERSConventions.IERS_2010, TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True),
                                    [CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()])



            # Atmospheric drag
            #from org.orekit.models.earth.atmosphere import NRLMSISE00
            #atmosphere = NRLMSISE00(msafe, sun, wgs84Ellipsoid)
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





    def ExBuildWalker(self, num_plane, num_sat, F, refSat):
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
                                          ta * ( np.pi / 180 ), PositionAngle.TRUE,
                                           refSat.getFrame(), refSat.getPVCoordinates().getDate(), 
                                           Constants.WGS84_EARTH_MU)

                orbits.append(newOrbit)
        return orbits

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


elevationMask = 10.0 # in degrees

# Scenario Start Time
utc = TimeScalesFactory.getUTC()
startTime = AbsoluteDate(2022, 10, 3, 8, 30, 0.0, utc)
endTime = AbsoluteDate(2022, 10, 3, 8, 40, 0.0, utc)
#referenceDate = AbsoluteDate.J2000_EPOCH  # reference date
#mjdUtcEpoch = AbsoluteDate(1858, 11, 17, 0, 0, 0.0, utc)


ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
eciFrame = FramesFactory.getGCRF()
moon = CelestialBodyFactory.getMoon()
sun = CelestialBodyFactory.getSun()
wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecefFrame)







initialRealOrbitSat1 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 500000.0, # Semimajor Axis (m)
                0.0 ,    # Eccentricity
                78.0*np.pi/180,  # Inclination (radians)
                0.0*np.pi/180,   # Perigee argument (radians)
                180.0*np.pi/180 ,   # Right ascension of ascending node (radians)
                25.0*np.pi/180,  # Anomaly (rad/s)
                PositionAngle.TRUE,  # Sets which type of anomaly we use
                eciFrame, # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                startTime,   # Sets the date of the orbital parameters
                Constants.WGS84_EARTH_MU)   # Sets the central attraction coefficient (m³/s²)
initialRealOrbitSat2 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 500000.0, # Semimajor Axis (m)
                0.0 ,    # Eccentricity
                78.0*np.pi/180,  # Inclination (radians)
                0.0*np.pi/180,   # Perigee argument (radians)
                207.0*np.pi/180 ,   # Right ascension of ascending node (radians)
                25.0*np.pi/180,  # Anomaly (rad/s)
                PositionAngle.TRUE,  # Sets which type of anomaly we use
                eciFrame, # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                startTime,   # Sets the date of the orbital parameters
                Constants.WGS84_EARTH_MU)   # Sets the central attraction coefficient (m³/s²)
initialEstimatedOrbitSat1 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 501234.0,
                0.0, 78.1*np.pi/180, 0.0*np.pi/180, 180.0*np.pi/180, 25.00*np.pi/180,
                PositionAngle.TRUE, eciFrame, startTime,   Constants.WGS84_EARTH_MU)
initialEstimatedOrbitSat2 = KeplerianOrbit(Constants.WGS84_EARTH_EQUATORIAL_RADIUS + 500221.0,
                0.0, 78.0*np.pi/180, 0.0*np.pi/180, 207.0*np.pi/180, 25.0*np.pi/180, 
                PositionAngle.TRUE, eciFrame, startTime, Constants.WGS84_EARTH_MU)


# Orbit propagator parameters
prop_min_step = 0.01 # s
prop_max_step = 2.0 # s
prop_position_error = 3.0 # m

sat1RealPropagator = COD_.ExGetPropagator("Get",initialRealOrbitSat1,100.0,
                                  eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                   prop_position_error, prop_position_error)
sat2RealPropagator = COD_.ExGetPropagator("Get",initialRealOrbitSat2,100.0,
                                  eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                   prop_position_error, prop_position_error)

sat1EstimatedPropagatorBuilder = COD_.ExGetPropagator("Builder",initialEstimatedOrbitSat1,100.0,
                                  eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                   prop_position_error, prop_position_error)

sat2EstimatedPropagatorBuilder = COD_.ExGetPropagator("Builder",initialEstimatedOrbitSat2,100.0,
                                  eciFrame, ecefFrame, "H", prop_min_step, prop_max_step,
                                   prop_position_error, prop_position_error)








estimator_convergence_thres = 1e-8
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

range_weight = 0.01 # Will be normalized later (i.e divided by the number of observations)
range_sigma = 0.01 # Estimated covariance of the range measurements, in meters

isTwoWay = True
currentDateTime = startTime

this_cons = COD_()
this_cons.addGroundStation(35.0, 51.0, 0.0, "Tehran")
this_cons.addGroundStation(30.0, 50.0, 0.0, "Boushehr")
this_cons.addGroundStation(33.0, 48.0, 1500.0, "Aleshtar")
this_cons.addGroundStation(25.0, 59.0, 0.0, "Bandar Abbas")
this_cons.addGroundStation(35.0, 58.0, 0.0, "Mashhad")

for i in range(2):
    meas = COD_.generate_measurements(sat1RealPropagator,this_cons.groundStationsList[i],str(i),0.0,5.0*np.pi/180.0,"RANGERATE",
                                    currentDateTime,1.0)
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
sat1EstimatedPropagator.propagate(startTime.shiftedBy(300.0), startTime.shiftedBy(2 * 3600.0 + 300))
sat1BoundedPropagator = sat1EphemerisGenerator.getGeneratedEphemeris()


currentDateTime = startTime.shiftedBy(300.0)
endTime = startTime.shiftedBy(2 * 3600.0 + 300.0)
position_norm_resi = pd.DataFrame(columns=['PositionNormError'])
while currentDateTime.compareTo(endTime) <= 0:
    x_resi = (sat1RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getX() - 
              sat1BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getX())
    y_resi = (sat1RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getY() - 
              sat1BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getY())
    z_resi = (sat1RealPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getZ() - 
              sat1BoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition().getZ())
    norm_resi = (x_resi**2 + y_resi**2 + z_resi**2)**0.5
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
    title = 'Position residuals',
    xaxis = dict(
        title = 'Datetime UTC'
    ),
    yaxis = dict(
        title = 'Position residual (m)'
    )
)

fig = dict(data=data, layout=layout)
#pio.write_image(fig, "file.png", height=1200, width=1600,scale=1)
pio.show(fig)

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
    norm_resi = (x_resi**2 + y_resi**2 + z_resi**2)**0.5
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
    title = 'Position residuals',
    xaxis = dict(
        title = 'Datetime UTC'
    ),
    yaxis = dict(
        title = 'Position residual (m)'
    )
)

fig = dict(data=data, layout=layout)

pio.show(fig)
# pio.write_image(fig, "test.svg", width=1.5*300, height=0.75*300, scale=1)
