% %
% % /* Copyright 2002-2019 CS Systèmes d'Information
% %  * Licensed to CS Systèmes d'Information (CS) under one or more
% %  * contributor license agreements.  See the NOTICE file distributed with
% %  * this work for additional information regarding copyright ownership.
% %  * CS licenses this file to You under the Apache License, Version 2.0
% %  * (the "License"); you may not use this file except in compliance with
% %  * the License.  You may obtain a copy of the License at
% %  *
% %  *   http://www.apache.org/licenses/LICENSE-2.0
% %  *
% %  * Unless required by applicable law or agreed to in writing, software
% %  * distributed under the License is distributed on an "AS IS" BASIS,
% %  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% %  * See the License for the specific language governing permissions and
% %  * limitations under the License.
% %  */
% %
% % Translated from SlaveMode.java to Matlab by Petrus Hyvönen 2011 as an
% % example of how to access orekit from matlab
% % the jars orekit-8.0.jar, hipparchus-core-1.0.jar, hipparchus-geometry-1.0.jar,
% % hipparchus-ode-1.0.jar, hipparchus-fitting-1.0.jar, hipparchus-optim-1.0.jar
% % and orekit-data.zip is in current matlab dir.
% 
% % These seems to work if pasted to prompt.
% javaaddpath 'D:\sabzy\MatlabToolBoxes\orekitInMatlab\orekit-12.0.jar'
% javaaddpath 'D:\sabzy\MatlabToolBoxes\orekitInMatlab\hipparchus-core-3.0.jar'
% javaaddpath 'D:\sabzy\MatlabToolBoxes\orekitInMatlab\hipparchus-geometry-3.0.jar'
% javaaddpath 'D:\sabzy\MatlabToolBoxes\orekitInMatlab\hipparchus-ode-3.0.jar'
% javaaddpath 'D:\sabzy\MatlabToolBoxes\orekitInMatlab\hipparchus-fitting-3.0.jar'
% javaaddpath 'D:\sabzy\MatlabToolBoxes\orekitInMatlab\hipparchus-optim-3.0.jar'
% 
% tic
% %% do the imports
% import org.orekit.utils.Constants
% import org.orekit.errors.OrekitException
% import org.orekit.frames.Frame
% import org.orekit.frames.FramesFactory
% import org.orekit.orbits.KeplerianOrbit
% import org.orekit.orbits.Orbit
% import org.orekit.orbits.PositionAngle
% import org.orekit.propagation.SpacecraftState
% import org.orekit.propagation.analytical.KeplerianPropagator
% import org.orekit.data.DataProvidersManager
% import org.orekit.data.ZipJarCrawler
% import org.orekit.time.AbsoluteDate
% import org.orekit.time.TimeScalesFactory
% import java.lang.Math;
% import java.lang.System;
% import java.util.Calendar;
% import java.util.Date;
% import java.util.GregorianCalendar;
% import java.util.TimeZone;
% import java.io.File;
% import org.orekit.models.earth.atmosphere.data.*;
% %% Configure Orekit. The file orekit-data.zip must be in current dir
% %Here is where things change in the setup
% DM=org.orekit.data.DataContext.getDefault().getDataProvidersManager(); %Works for v11.1
% %DM      = org.orekit.data.DataProvidersManager(); %throws error foe v11.1
% crawler = org.orekit.data.DirectoryCrawler(File('D:\sabzy\MatlabToolBoxes\orekitInMatlab\'));
% %crawler = org.orekit.data.DirectoryCrawler(File([cd ,'\OrekitData\']));
% 
% DM.clearProviders()
% DM.addProvider(crawler)
% %% Initial orbit parameters
% a = 6378137+700000;    % semi major axis in meters
% e = 0.00072831215;  % eccentricity
% i = (58/180)*pi; % inclination
% omega = (180)/180*pi; % perigee argument
% raan = (261)/180*pi; %right ascension of ascending node
% lM = 0.0; % mean anomaly
% a = MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES;
% b = javaMethod('valueOf','org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation$StrengthLevel', 'AVERAGE');
% %% Set inertial frame
% inertialFrame = FramesFactory.getEME2000();
% 
% %% Initial date in UTC time scale
% utc = TimeScalesFactory.getUTC();
% initialDate = AbsoluteDate(2004, 01, 01, 23, 30, 00.000, utc);
% 
% %% Setup orbit propagator
% %gravitation coefficient
% mu =  3.986004415e+14;
% 
% %Orbit construction as Keplerian
% initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN, inertialFrame, initialDate, mu);
% 
% %Simple extrapolation with Keplerian motion
% kepler = KeplerianPropagator(initialOrbit);
% 
% %Set the propagator to slave mode (could be omitted as it is the default mode)
% %kepler.setSlaveMode()
% 
% %% Setup propagation time
% %Overall duration in seconds for extrapolation
% duration = 180*60.0;
% 
% %Stop date
% finalDate =  AbsoluteDate(initialDate, duration, utc);
% 
% %Step duration in seconds
% stepT = 30.0;
% 
% %% Perform propagation
% %Extrapolation loop
% cpt = 1;
% extrapDate = initialDate;
% while extrapDate.compareTo(finalDate) <= 0
%     currentState = kepler.propagate(extrapDate);
%     %fprintf('step %d: time %s %s\n', cpt, char(currentState.getDate()), char(currentState.getOrbit()))
%     coord=currentState.getPVCoordinates.getPosition;
%     P(:,cpt)=[coord.getX coord.getY coord.getZ]';
%     extrapDate = AbsoluteDate(extrapDate, stepT, utc);
%     cpt=cpt+1;
% end
% %%
% toc
% [X,Y,Z] = sphere(50);
% figure;
% plot3(P(1,:),P(2,:),P(3,:),'LineWidth',2);
% hold on
% surf(X.*6378137,Y.*6378137,Z.*6378137)
% axis equal
%
