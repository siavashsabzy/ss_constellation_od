import orekit
vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
setup_orekit_curdir()

from orekit import JArray, JArray_double

# Measurements Builders
from org.orekit.estimation.measurements.generation import Generator
from org.orekit.estimation.measurements.generation import RangeBuilder
from org.orekit.estimation.measurements.generation import RangeRateBuilder
from org.orekit.estimation.measurements.generation import AngularAzElBuilder

# Modifiers - Measurements Perturbations
from org.orekit.estimation.measurements.modifiers import RangeTroposphericDelayModifier
from org.orekit.estimation.measurements.modifiers import RangeRateTroposphericDelayModifier

from org.orekit.estimation.measurements.modifiers import RangeIonosphericDelayModifier
from org.orekit.estimation.measurements.modifiers import RangeRateIonosphericDelayModifier

from org.orekit.estimation.measurements.modifiers import AngularTroposphericDelayModifier
from org.orekit.estimation.measurements.modifiers import AngularIonosphericDelayModifier

from org.orekit.estimation.measurements.modifiers import TurnAroundRangeIonosphericDelayModifier
from org.orekit.estimation.measurements.modifiers import TurnAroundRangeTroposphericDelayModifier

from org.orekit.estimation.measurements import Range
from org.orekit.estimation.measurements import RangeRate
from org.orekit.estimation.measurements import AngularAzEl
from org.orekit.estimation.measurements import GroundStation
from org.orekit.estimation.measurements import ObservableSatellite

from org.orekit.utils import Constants
from org.orekit.bodies import GeodeticPoint
from org.orekit.frames import TopocentricFrame
from org.orekit.estimation.measurements import GroundStation
# Atmospheric Models
from org.orekit.models.earth import EarthITU453AtmosphereRefraction
from org.orekit.models.earth.weather import GlobalPressureTemperature2Model
# Tropospheric Models
from org.orekit.models.earth.troposphere import SaastamoinenModel
from org.orekit.models.earth.troposphere import MariniMurrayModel # Mainly for laser Ranging
# Ionospheric Models
from org.orekit.models.earth.ionosphere import EstimatedIonosphericModel
from org.orekit.models.earth.ionosphere import GlobalIonosphereMapModel
from org.orekit.models.earth.ionosphere import KlobucharIonoModel
from org.orekit.models.earth.ionosphere import NeQuickModel

from org.orekit.data import DataProvidersManager, ZipJarCrawler
from org.orekit.time import TimeScalesFactory, AbsoluteDate
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation import SpacecraftState
from org.orekit.orbits import OrbitType, PositionAngle
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.hipparchus.ode.nonstiff import ClassicalRungeKuttaIntegrator

from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.gravity import ThirdBodyAttraction
from org.orekit.forces.gravity import Relativity

from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
from org.orekit.forces.radiation import SolarRadiationPressure

from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.utils import Constants
from org.orekit.utils import IERSConventions

from org.orekit.utils import Constants, PVCoordinates, PVCoordinatesProvider, IERSConventions
#from org.orekit.estimation import Context
from org.orekit.estimation.leastsquares import BatchLSEstimator, PythonBatchLSObserver
from org.orekit.estimation.measurements import GroundStation, Range, RangeRate, AngularAzEl, PV, ObservedMeasurement, ObservableSatellite
from org.orekit.estimation.measurements.modifiers import Bias, OutlierFilter
from org.orekit.estimation.measurements.generation import EventBasedScheduler, SignSemantic
from org.hipparchus.optim.nonlinear.vector.leastsquares import LevenbergMarquardtOptimizer, GaussNewtonOptimizer
from org.hipparchus.linear import QRDecomposer, Array2DRowRealMatrix
from org.orekit.bodies import GeodeticPoint, OneAxisEllipsoid, CelestialBodyFactory
from org.orekit.models.earth.displacement import StationDisplacement
from org.orekit.time import TimeScalesFactory, AbsoluteDate, DateTimeComponents, FixedStepSelector
import org.orekit.time as oktime
from org.orekit.orbits import CartesianOrbit, KeplerianOrbit, EquinoctialOrbit, PositionAngle, OrbitType
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel, OceanTides, Relativity, SolidTides, ThirdBodyAttraction, NewtonianAttraction
from org.orekit.forces.gravity.potential import GravityFieldFactory, NormalizedSphericalHarmonicsProvider
from org.orekit.propagation import Propagator, SpacecraftState
from org.orekit.propagation.conversion import NumericalPropagatorBuilder, DormandPrince853IntegratorBuilder
from org.orekit.propagation.numerical import NumericalPropagator, PartialDerivativesEquations
from org.orekit.propagation.events import ElevationDetector
from org.orekit.propagation.events.handlers import ContinueOnEvent

from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator

import sys
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from math import degrees, radians

class Orbit:
	def __init__(self, apogee, perigee, inclination, aop, raan, lv):
		self.apogee         = apogee*1000
		self.perigee        = perigee*1000
		self.inclination    = radians(inclination)
		self.aop            = radians(aop)
		self.raan           = radians(raan)
		self.lv             = radians(lv)

	def sma(self):
		return (self.perigee + self.apogee + 2*Constants.WGS84_EARTH_EQUATORIAL_RADIUS)/2.0

	def eccentricity(self):
		return abs((self.apogee - self.perigee)/(self.apogee + self.perigee + 2*Constants.WGS84_EARTH_EQUATORIAL_RADIUS))

class Spacecraft:
	def __init__(self, name, mass, cross_section, cd, cr, identification):
		self.name           = name
		self.mass           = mass
		self. cross_section = cross_section
		self.cd             = cd
		self.cr             = cr
		self.identification = identification

class GroundSegment:
	def __init__(self, name, latitude, longitude, altitude, elevationAngle, temperature, humidity, pressure, identification):
		self.name           = name
		self.latitude       = latitude
		self.longitude      = longitude
		self.altitude       = altitude
		self.elevationAngle = radians(elevationAngle)
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

def numericalPropagator(initialDate, sma, ecc, inc, aop, raan, lv, mass, cross_section,
	cr     = 1.8,
	degree = 70,
	order  = 70): 

	inertialFrame = FramesFactory.getEME2000()
	terrestrialFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
	earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, terrestrialFrame)

	initialOrbit = KeplerianOrbit(sma, ecc, inc, aop, raan, lv, PositionAngle.TRUE, inertialFrame, initialDate, Constants.WGS84_EARTH_MU)
	initialState = SpacecraftState(initialOrbit, mass)
	
	# Setup numerical propagator
	minStep  = 1e-3
	maxStep  = 1e+3
	initStep = 60.0

	positionTolerance = 1.0
	orbitType = OrbitType.KEPLERIAN
	tolerance = NumericalPropagator.tolerances(positionTolerance, initialOrbit, orbitType)
	
	# Integrator
	integrator = DormandPrince853Integrator(minStep, maxStep,
		JArray_double.cast_(tolerance[0]),
		JArray_double.cast_(tolerance[1]))
	integrator.setInitialStepSize(initStep)

	propagator = NumericalPropagator(integrator)
	propagator.setOrbitType(orbitType)
	
	# GravityModel
	gravityModel = GravityFieldFactory.getNormalizedProvider(degree, order)
	propagator.addForceModel(HolmesFeatherstoneAttractionModel(earth.getBodyFrame(), gravityModel))

	# Third-Body Forces - Sun and Moon
	sun3rdBodyAttraction  = ThirdBodyAttraction(CelestialBodyFactory.getSun())
	moon3rdBodyAttraction = ThirdBodyAttraction(CelestialBodyFactory.getMoon())

	# Solar Radiation Pressure
	isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(cross_section, cr)
	solarRadiationPressure = SolarRadiationPressure(CelestialBodyFactory.getSun(),
		Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
		isotropicRadiationSingleCoeff)

	# Relativity Correction
	relativityCorrection = Relativity(Constants.WGS84_EARTH_MU)

	# Add Perturbations
	propagator.addForceModel(sun3rdBodyAttraction)
	propagator.addForceModel(moon3rdBodyAttraction)
	propagator.addForceModel(solarRadiationPressure)
	propagator.addForceModel(relativityCorrection)

	propagator.setInitialState(initialState)

	return propagator

class StationError:
	def __init__(self,azelSigma, rangeSigma, rangeRateSigma):
		self.azelSigma      = azelSigma
		self.rangeSigma     = rangeSigma
		self.rangeRateSigma = rangeRateSigma

def orekitSpaceInventorStations(stations):
	groundstations = []
	groundstationNames = []

	for station in stations:
		groundstations.append(station)
		groundstationNames.append(station.name)

	return groundstations, groundstationNames

def noise_generator(seed, sigma, small = 1e-10):
    from org.hipparchus.random import Well19937a, GaussianRandomGenerator, CorrelatedRandomVectorGenerator
    from org.hipparchus.linear import MatrixUtils
    
    random_generator    = Well19937a(int(seed))
    gaussian_generator  = GaussianRandomGenerator(random_generator)
    #covariance = MatrixUtils.createRealDiagonalMatrix(float(sigma * sigma))
    covariance = MatrixUtils.createRealDiagonalMatrix([float(sigma * sigma),float(sigma * sigma)])
    
    return CorrelatedRandomVectorGenerator(covariance, float(small), gaussian_generator)

def generate_measurements(propagator, station, station_name, altitude, elevAngle, meas_type, t0, duration,
	sigma          = None,
	base_weight    = 1.0,
	two_way        = True,
	withRefraction = True,
	step           = 3600.0,
	seed           = 0):

	tf = t0.shiftedBy(3600.0*duration)

	noise_source = noise_generator(seed, sigma = sigma)

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
		if two_way:
			way_str = "TWOWAY"
		else:
			way_str = "ONEWAY"
		print("{}   {} {}       {}        {}".format(meas.getDate().toString(), way_str, meas_type, station_name, meas.getObservedValue()[0]))


if __name__ == '__main__':

	utc = TimeScalesFactory.getUTC()

	# Frame and Earth definitions
	itrf93 = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
	earth  = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, itrf93)

	# Beginning Of Life (BOL)
	date = datetime(2022, 2, 1, 10, 30, 0)  # (YYYY,M,D,H,M,S)
	# Convert to OREKIT Date
	epochDate   = AbsoluteDate(date.year,date.month,date.day,date.hour,date.minute,float(date.second),utc)
	initialDate = epochDate

	# Initial State
	gravitySpacecraft = Spacecraft("GravitySpace", 21.86, 0.29731552, 2.2, 1.8, 11111)
	gravityOrbit      = Orbit(814, 786, 98.55, 90.0, 5.1917, 0.0567634)

	# Setup numerical propagator
	gravityProp = numericalPropagator(initialDate, gravityOrbit.sma(), gravityOrbit.eccentricity(),
		gravityOrbit.inclination, gravityOrbit.aop, gravityOrbit.raan, gravityOrbit.lv,
		gravitySpacecraft.mass, gravitySpacecraft.cross_section, gravitySpacecraft.cr, 70, 70)

	# Define SSC Ground Stations
	alaska      = GroundSegment("Alaska", +64.8000, -147.5000, 6.1905, 5.0, 295.1, 0.55, 1013.25, 22222)
	puntaArenas = GroundSegment("PuntaArenas", -52.9400, -70.8600, 0.0340, 5.0, 295.1, 0.55, 1013.25, 33333)
	esrange     = GroundSegment("Esrange", +67.8800, +21.0500, 0.5300, 5.0, 295.1, 0.55, 1013.25, 44444)
	dongara     = GroundSegment("Dongara", -29.0500, +115.3500, 0.0340, 5.0, 295.1, 0.55, 1013.25, 55555)

	stations    = (alaska, puntaArenas, esrange, dongara)
	#stations, stationNames = orekitSpaceInventorStations(stations)

	# Measurements Errors - Applicable to all ground stations from SSC
	sscSigma = StationError(0.080, 15.000, 0.040) # [deg, deg, m, m.s-1]

	# Ranging Data Generation
	base_weight = 1.0
	seed = int(0) # int(sys.argv[0])
	step = 3.0 # 60*float(sys.argv[2])
	#if sys.argv[3] == '1':
	rangeTwoWay = True
	#else:
	#	rangeTwoWay = False
	#if sys.argv[4] == '1':
	rateTwoWay = True
	#else:
	#	rateTwoWay = False

	for station in stations:

		#generate_measurements(gravityProp, station.groundStation(earth), station.name, station.altitude, station.elevationAngle,
		#	'RANGE', initialDate, 24.0*10, sigma=sscSigma.rangeSigma, two_way=rangeTwoWay, withRefraction=True, step=step, seed=seed)
		generate_measurements(gravityProp, station.groundStation(earth), station.name, station.altitude, station.elevationAngle,
			'RANGERATE', initialDate, 24.0*10, sigma=sscSigma.rangeRateSigma, two_way=rangeTwoWay, withRefraction=True, step=step, seed=seed)
		#generate_measurements(gravityProp, station.groundStation(earth), station.name, station.altitude, station.elevationAngle,
		#	'AZEL', initialDate, 24.0, sigma=sscSigma.azelSigma, two_way=rangeTwoWay, withRefraction=True, step=step, seed=seed)