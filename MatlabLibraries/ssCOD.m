classdef ssCOD
    % Copyright © 2002-2024 CS GROUP
    % @author: Siavash Sabzy (signature: ss...)
    % https://github.com/siavashsabzy

    properties
        startTime
        duration
        durationAfterEstimation
        endEstimation
        endTime
        plotStepTime
        RAD2DEG
        DEG2RAD
        nominalAltitudes
        showActiveSatellites
        rangeMaxError
        rangeSigma
        rangeWeight
        stationSigma
        stationBaseWeight
        stationisTwoWay
        stationRefractio
        stationSampleTime
        stationStd
        noiseSeed
        isISL
        islSigma
        islBaseWeight
        isTwoWay
        islSampleTime
        islScheduler
        islStd
        gpsSampleTime
        gpsSigma
        gpsBaseWeight
        gpsPositionStd
        gpsVelocitystd
        GPSMaxBurstSize
        GPSHighRateStep
        GPSBurstPeriod
        groundStationsList
        numberOfGroundStation
        spacecraftsList
        propagatorsList
        referenceOrbit
        ecefFrame
        eciFrame
        wgs84Ellipsoid
        moon
        sun
        propagationMinTimeStep
        propagationMaxTimeStep
        vectorAbsoluteTolerance
        vectorRelativeTolerance
        BLSestimatorConvergenceThreshold
        BLSestimatorMaxIterations
        BLSestimatorMaxEvaluations
        elevationMask
        mass
        islWaveLength
    end


    methods (Static)

        function addOrekitLibraries()
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\orekit-12.0.jar'
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-clustering-3.0.jar'
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-core-3.0.jar'
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-fft-3.0.jar'
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-filtering-3.0.jar';
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-fitting-3.0.jar'
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-geometry-3.0.jar';
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-migration-3.0.jar';
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-ode-3.0.jar';
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-optim-3.0.jar';
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-samples-3.0.jar';
            javaaddpath 'C:\sabzy\ss_constellation_od\MatlabLibraries\hipparchus-stat-3.0.jar';
        end

        function lla = ecefToLla(X)

            import org.orekit.utils.*;
            R_e = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
            Earth_E2 = 8.181919084261345e-2;
            r_I = X(1);
            r_J = X(2);
            r_K = X(3);
            if r_I == 0 && r_J == 0 % you are at one of the poles

                long    = 0;

                % polar radius of the Earth:
                r_p = R_e*sqrt(1-Earth_E2);

                if r_K == 0 % at center
                    h_ellp = -R_e;
                    lat = 0;
                elseif r_K > 0 % in Northern hemisphere
                    h_ellp = r_K - r_p;
                    lat = pi/2;
                else % r_K <0, in southern hemisphere
                    h_ellp = -(r_K + r_p);
                    lat = -pi/2;
                end

            else

                r_delta = sqrt((r_I^2)+(r_J^2));

                % compute alpha
                s_alp = r_J/r_delta;
                c_alp = r_I/r_delta;

                alpha = atan2(s_alp,c_alp);

                % compute delta
                delta = atan(r_K/r_delta);

                % compute longitude
                long = alpha;

                % compute geodetic latitude:
                %
                phi     = delta;   % initial guess
                phi_old = 1e5;     % dummy value
                tol     = 1e-10;    % tolerence to be used in the iteration

                % iterate
                while abs(phi_old-phi)>tol

                    temp = phi;
                    C_E     = R_e/sqrt(1-Earth_E2*sin(temp)^2);
                    phi     = atan((r_K+C_E*Earth_E2*sin(temp))/r_delta);
                    phi_old = temp;

                end

                lat = phi;

                % Using the values of phi_gd and C_E determined above, compute h_ellp
                h_ellp = (r_delta/cos(phi)-C_E);

            end
            lla = [lat,long,h_ellp];
        end

        function Rcv = llaToEcef(Stat)
            latit  =  Stat(1);   % ground station latitude
            longi  =  Stat(2);   % ground station longitude
            altit  =  Stat(3);   % ground station height
            e2     =  0.00669437999014132;   % Square of eccentricity
            CosLat =  cosd( latit );         % (Co)sine of geodetic latitude
            SinLat =  sind( latit );
            N      =  6378137.0 / sqrt( 1.0 - e2 * SinLat * SinLat );
            Rcv(1) =  ( N + altit) * CosLat * cosd( longi );
            Rcv(2) =  ( N + altit) * CosLat * sind( longi );
            Rcv(3) =  ( ( 1.0 - e2 ) * N + altit ) * SinLat;
        end

        function logarithms = thisLogarithm(values, bases)
            logarithms = log(values) ./ log(bases);
        end

        function subTriangles = divideTriangles(triangles, numDimensions, ...
                numVerticesInATriangle, numNewTrianglesPerTriangle)
            oldAs = triangles(:, :, 1);
            oldBs = triangles(:, :, 2);
            oldCs = triangles(:, :, 3);
            Ps = (oldAs + oldBs)/ 2;
            Qs = (oldBs + oldCs)/ 2;
            Rs = (oldCs + oldAs)/ 2;
            % Find the midpoints of each side.

            scaling = 1 / norm(Ps(1, :));
            unitPs = scaling * Ps;
            unitQs = scaling * Qs;
            unitRs = scaling * Rs;
            % Normalize midpoints onto the surface of a unit sphere.

            newAs = [oldAs; unitPs; unitRs; unitQs];
            newBs = [unitPs; oldBs; unitQs; unitRs];
            newCs = [unitRs; unitQs; oldCs; unitPs];
            % Find the sub-triangles' vertices. Ensure that every point gets used as an
            % A, B, and C in at least one triangle.

            subTriangles = reshape([newAs, newBs, newCs], ...
                numNewTrianglesPerTriangle * uint32(size(triangles, 1)), numDimensions, ...
                numVerticesInATriangle);
            % Put the sub-triangles' vertices in an array of the same form as the
            % original triangles.

        end

        function [latsInDegrees, longsInDegrees] = cartesianToSpherical(Xs, Ys, Zs)
            latsInDegrees = asind(Zs ./ realsqrt((Xs .^ 2) + (Ys .^ 2) + (Zs .^ 2)));
            centers = ~isfinite(latsInDegrees);
            latsInDegrees(centers) = NaN;
            longsInDegrees = (180 / pi())*(atan2(Ys, Xs));
            longsInDegrees(centers) = NaN;
        end

        function distinctFromPrevious = gridDistinctFromPrevious(values)
            shiftedVals = circshift(values, 1);
            if (isinteger(values))
                equalsPrevious = (values == shiftedVals);
            else
                TOLERANCE_FACTOR = 32;
                maxMagnitudes = max(abs(values), abs(shiftedVals));
                tolerance = TOLERANCE_FACTOR * arrayfun(@eps, maxMagnitudes);
                equalsPrevious = (abs(values - shiftedVals) <= tolerance);
            end
            equalsPrevious(1) = false;
            distinctFromPrevious = uint32(find(~equalsPrevious));
        end

        function [GDOP, PDOP, TDOP, HDOP, VDOP] = getSiteDops(R_site, ECEF2ENU, R_sat, el_mask)
            GDOP = NaN;PDOP = NaN;TDOP = NaN;HDOP = NaN;VDOP = NaN;
            R_e = 6378137;
            sin_el_mask = sin(el_mask*pi/180);
            radius_sat_squared = sum( R_sat.^2, 2 );
            rho_mask = ...
                -R_e * sin_el_mask +...
                (R_e ^ 2 * (sin_el_mask ^ 2 - 1) + radius_sat_squared).^0.5;

            test = R_sat*R_site';
            metric = ( R_e^2 + radius_sat_squared - rho_mask.^2 )./2;
            sv_in_view = test>=metric;
            num_sv = sum(sv_in_view);
            if num_sv >= 4
                R_sat_hat = R_sat(sv_in_view,:);
                LOS = R_sat_hat - repmat(R_site, num_sv, 1);
                rho = sum( LOS.^2, 2).^0.5;
                RHO = repmat(rho, 1, 3);
                G = [LOS./RHO, ones(num_sv,1)];
                R_l_tilde = eye(4);
                R_l_tilde(1:3, 1:3) = ECEF2ENU;
                H = inv(R_l_tilde*(G'*G)*R_l_tilde');
                DOPs = diag(H);
                GDOP = sqrt( sum(DOPs) );
                PDOP = sqrt( sum(DOPs(1:3)) );
                TDOP = sqrt( DOPs(4) );
                HDOP = sqrt( sum(DOPs(1:2)) );
                VDOP = sqrt( DOPs(3) );
            end


        end

        function [coe, undefined, orbit_type] = ECI2COE(X_ECI,V_ECI)
            import org.orekit.utils.*;
            mu = Constants.WGS84_EARTH_MU;
            r = norm(X_ECI);
            v = norm(V_ECI);
            h_vec = cross(X_ECI,V_ECI);
            h     = norm(h_vec);
            n_vec = cross([0 0 1]',h_vec);
            n     = norm(n_vec);
            energy = (v^2)/2 - mu/r;
            e_vec = ((v^2 - mu/r)*X_ECI - dot(X_ECI,V_ECI)*V_ECI ) / mu;
            e = norm(e_vec);
            tol = 1e-15;
            if abs(e-1)<tol
                e = 1;
            end
            if abs(e)<tol
                e = 0;
            end
            if e ~= 1
                a = -mu/2/energy;
                p = a*(1-e^2);
            else % parabolic case
                a = Inf;
                p = h^2/mu;
            end
            i = acos(h_vec(3)/h); % [rad]
            RAAN = acos( n_vec(1)/n );
            if n_vec(2)<0
                RAAN = 2*pi-RAAN;
            end
            w = real(acos( dot(n_vec,e_vec) / e / n )); % [rad]
            if e_vec(3)<0
                w = 2*pi-w;
            end
            f = real(acos( dot(e_vec,X_ECI) / r / e )); % [rad]
            if dot(X_ECI,V_ECI)<0
                f = 2*pi-f;
            end
            cosf = cos(f);
            cosE = (e+cosf)/(1+e*cosf);
            sinE = sin(f)*sqrt(1-e^2)/(1+e*cosf);
            E = atan2(sinE,cosE);
            M = E - e*sin(E);

            orbit_case = 1; % Default.

            if (e<1 && e~=0) && (i == 0 || i == pi)
                orbit_case = 2;
            end

            if e == 0 && (i ~= 0 || i ~= pi)
                orbit_case = 3;
            end

            if e == 0 && (i == 0 || i == pi)
                orbit_case = 4;
            end

            switch orbit_case
                case 1
                    coe.p = p;
                    coe.a = a;
                    coe.e = e;
                    coe.i = i*180/pi;
                    coe.omega = w*180/pi;
                    coe.RAAN = RAAN*180/pi;
                    coe.f = f*180/pi;
                    coe.M = M*180/pi;
                    undefined.p = 0;
                    undefined.a = 0;
                    undefined.e = 0;
                    undefined.i = 0;
                    undefined.omega = 0;
                    undefined.RAAN = 0;
                    undefined.f = 0;
                    undefined.M = 0;
                    orbit_type = 'elliptical inclined';

                case 2
                    w_true = acos(e_vec(1)/e); % [rad]
                    if e_vec(2)<0
                        w_true =  2*pi-w_true;
                    end
                    coe.p = p;
                    coe.a = a;
                    coe.e = e;
                    coe.i = i*180/pi;
                    coe.omega = w_true*180/pi;
                    coe.RAAN = 0;
                    coe.f = f*180/pi;
                    coe.M = M*180/pi;
                    undefined.p = 0;
                    undefined.a = 0;
                    undefined.e = 0;
                    undefined.i = 0;
                    undefined.omega = 1;
                    undefined.RAAN = 1;
                    undefined.f = 0;
                    undefined.M = 0;
                    orbit_type = 'elliptical equitorial';

                case 3
                    u = acos( dot(n_vec,X_ECI) / n / r );
                    if X_ECI(3)<0
                        u = 2*pi-u;
                    end
                    coe.p = p;
                    coe.a = a;
                    coe.e = e;
                    coe.i = i*180/pi;
                    coe.omega = 0;
                    coe.RAAN = RAAN*180/pi;
                    coe.f = u*180/pi;
                    coe.M = coe.f;
                    undefined.p = 0;
                    undefined.a = 0;
                    undefined.e = 0;
                    undefined.i = 0;
                    undefined.omega = 1;
                    undefined.RAAN = 0;
                    undefined.f = 1;
                    undefined.M = 1;
                    orbit_type = 'circular inclined';

                case 4
                    lambda_true = acos(X_ECI(1) / r);
                    if X_ECI(2)<0
                        lambda_true = 2*pi-lambda_true;
                    end
                    coe.p = p;
                    coe.a = a;
                    coe.e = e;
                    coe.i = i*180/pi;
                    coe.omega = 0;
                    coe.RAAN = 0;
                    coe.f = lambda_true*180/pi;
                    coe.M = coe.f;
                    undefined.p = 0;
                    undefined.a = 0;
                    undefined.e = 0;
                    undefined.i = 0;
                    undefined.omega = 1;
                    undefined.RAAN = 1;
                    undefined.f = 1;
                    undefined.M = 1;
                    orbit_type = 'circular equitorial';
            end
        end

        function [coeff_A2, coeff_C2, coeff_R2, coeff_T2, coeff_RT, theta] = ...
                analytic_URE_eqn( altitude )
            import org.orekit.utils.*;
            R_e = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
            digits(5);
            % Compute the orbital radius.
            R_s = ( R_e + altitude ); % []
            % Compute theta.
            theta = asin(R_e / R_s); % [rad]
            % disp(['theta = ', num2str(theta*180/pi),' [deg]'])
            % Define symbolic variables.
            alpha = sym('alpha');
            beta = sym('beta');
            T = sym('T');
            C = sym('C');
            A = sym('A');
            R = sym('R');
            % Add constraints on these symbolic variables.
            assume(alpha,'real')
            assume(beta,'real')
            assume(T,'real')
            assume(C,'real')
            assume(A,'real')
            assume(R,'real')
            assume(R,'positive')
            integrand2 = -(sin(alpha)*(23450783740712698330090463749682*R^2*cos(alpha)^2 -...
                11725391870356349165045231874841*C^2*cos(alpha)^2 - 11725391870356349165045231874841*A^2*cos(alpha)^2 ...
                + 7353490718638634059169792*R^2*altitude + 7353490718638634059169792*T^2*altitude +...
                11725391870356349165045231874841*A^2 + 11725391870356349165045231874841*C^2 ...
                + 23450783740712698330090463749682*R^2 + 46901567481425396940485430273586*T^2 + ...
                576460752303423488*R^2*altitude^2 + 576460752303423488*T^2*altitude^2 - ...
                46901567481425396660180927499364*R^2*cos(alpha) - 46901567481425396660180927499364*T^2*cos(alpha) ...
                - 7353490718638634059169792*R^2*altitude*cos(alpha) - 7353490718638634059169792*T^2*altitude*cos(alpha) ...
                - 7353490718638634059169792*R*T*((3424235954246779*altitude)/268435456 - ...
                (11725391870356349165045231874841*cos(alpha))/144115188075855872 - ...
                (3424235954246779*altitude*cos(alpha))/268435456 + altitude^2 + ...
                23450783740712698470242715136793/288230376151711744)^(1/2) - ...
                1152921504606846976*R*T*altitude*((3424235954246779*altitude)/268435456 - ...
                (11725391870356349165045231874841*cos(alpha))/144115188075855872 - ...
                (3424235954246779*altitude*cos(alpha))/268435456 + altitude^2 + ...
                23450783740712698470242715136793/288230376151711744)^(1/2) + ...
                7353490718638634059169792*R*T*cos(alpha)*((3424235954246779*altitude)/268435456 - ...
                (11725391870356349165045231874841*cos(alpha))/144115188075855872 - ...
                (3424235954246779*altitude*cos(alpha))/268435456 + altitude^2 + ...
                23450783740712698470242715136793/288230376151711744)^(1/2)))/(2*(3424235954246779/(536870912*(altitude ...
                + 3424235954246779/536870912)) - 1)*(3676745359319317029584896*altitude - ...
                23450783740712698330090463749682*cos(alpha) - 3676745359319317029584896*altitude*cos(alpha) + ...
                288230376151711744*altitude^2 + 23450783740712698470242715136793));
            % Integrate with respect to alpha.
            F = int(integrand2, alpha, 0, pi/2 - theta);

            %Display the equation.
            fprintf('%s','   ')
            equee = vpa(F);
            fprintf('%s','URE equation is:  ')
            fprintf('%s\n',equee)
            fprintf('%s\n','   ')
            fprintf('%s','              ')
            fprintf('%s\n','       ***************[Based on GPS SPS]*****************')

            % Get the coefficients of the polynomial.
            [coeff, ~] = coeffs(vpa(F),A);
            coeff_A2 = double( coeff(1) );

            [coeff, ~] = coeffs(vpa(F),C);
            coeff_C2 = double( coeff(1) );

            [coeff, ~] = coeffs(vpa(F),R);
            coeff_R2 = double( coeff(1) );

            [coeff, ~] = coeffs(vpa(F),T);
            coeff_T2 = double( coeff(1) );

            [coeff, ~] = coeffs(coeff(2),R);
            coeff_RT = double( coeff(1) );
        end

        function E_a = Keplers_Eqn(M, e)
            % Select initial guess.
            if ( (-pi < M) && (M < 0) ) || (M > pi)
                E_a = M-e;
            else
                E_a = M+e;
            end
            % Define tolerance.
            tol  = 1e-12;
            test = 999; % Dummy variable.
            % Implement Newton's method.
            while test > tol
                E_new = E_a + (M-E_a+e*sin(E_a))/(1-e*cos(E_a));
                test = abs(E_new - E_a);
                E_a = E_new;
            end
        end

        function [X_ECI,V_ECI] = ECEF2ECI(X_ECEF,V_ECEF,GMST)
            import org.orekit.utils.*;
            omega_e = Constants.WGS84_EARTH_ANGULAR_VELOCITY;
            ECI_C_ECEF = [cos(GMST) -sin(GMST) 0;sin(GMST) cos(GMST) 0;0 0 1] ;
            X_ECI = ECI_C_ECEF*X_ECEF;
            V_ECI = ECI_C_ECEF*V_ECEF + cross([0 0 omega_e]',X_ECI);
        end

        function ECI_2_RIC_Mat = ECI2RIC( x_eci, v_eci )
            % Radial.
            r_hat = x_eci / norm(x_eci);
            % Cross Track.
            c_hat = cross( x_eci, v_eci );
            c_hat = c_hat / norm( c_hat );
            % In / Along Track.
            i_hat = cross( c_hat, r_hat );
            % Return the transformation matrix.
            ECI_2_RIC_Mat = [ r_hat'; i_hat'; c_hat'];
        end

        function [N1,N2] = trilateration(P,S,W)
            [mp,np] = size(P);
            ns = length(S);
            if (ns~=np)
                error('Number of reference points and distances are different');
            end
            A=[]; b=[];
            for i1=1:np
                x = P(1,i1); y = P(2,i1); z = P(3,i1);
                s = S(i1);
                A = [A ; 1 -2*x  -2*y  -2*z];
                b= [b ; s^2-x^2-y^2-z^2 ];
            end
            if (np==3)
                warning off;
                Xp= A\b;  % Gaussian elimination
                % or Xp=pinv(A)*b;
                % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
                % depend only on the reference points
                % it could be computed only once
                xp = Xp(2:4,:);
                Z = null(A,'r');
                z = Z(2:4,:);
                if rank (A)==3
                    %Polynom coeff.
                    a2 = z(1)^2 + z(2)^2 + z(3)^2 ;
                    a1 = 2*(z(1)*xp(1) + z(2)*xp(2) + z(3)*xp(3))-Z(1);
                    a0 = xp(1)^2 +  xp(2)^2+  xp(3)^2-Xp(1);
                    p = [a2 a1 a0];
                    t = roots(p);

                    %Solutions
                    N1 = Xp + t(1)*Z;
                    N2 = Xp + t(2)*Z;
                end
            end
            if  (np>3)
                %Particular solution

                if W~=diag(ones(1,length(W)))
                    C = W'*W;
                    Xpdw =inv(A'*C*A)*A'*C*b; % Solution with Weights Matrix
                else
                    Xpdw=pinv(A)*b; % Solution without Weights Matrix
                end

                % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
                % depend only on the reference points
                % it could be computed only once
                N1 = Xpdw;
                N2 = N1;
            end
        end

        function outputTime = getUnixTime(orekitTime)
            import org.oreki.time.*;
            outputTime = 0;
            dateTimeComps = orekitTime.getComponents(0);
            thisTime = dateTimeComps.getTime();
            thisDate = dateTimeComps.getDate();
            thisYear = thisDate.getYear();
            thisMonth = thisDate.getMonth();
            thisDay = thisDate.getDay();
            thisHour = thisTime.getHour();
            thisMinute = thisTime.getMinute();
            thisSecond = thisTime.getSecond();
            outputTime = datenum([thisYear, thisMonth, thisDay, thisHour, thisMinute, thisSecond]);
        end
    end


    methods

        function obj = ssCOD(inputStructure)
            import org.orekit.time.*;
            import java.io.File;
            import org.orekit.data.*;
            import org.orekit.frames.*;
            import org.orekit.models.earth.*;
            import org.orekit.bodies.*;
            import org.orekit.utils.*;
            import org.orekit.orbits.*;
            import org.orekit.gnss.*;
            DM = org.orekit.data.DataContext.getDefault().getDataProvidersManager();
            crawler = org.orekit.data.DirectoryCrawler(File('C:\sabzy\sabzyJava\OrekitLib\'));
            DM.clearProviders()
            DM.addProvider(crawler)
            obj.ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, ...
                IERSConventions.IERS_2010, true);
            obj.eciFrame = FramesFactory.getGCRF();
            obj.RAD2DEG = (180 / pi);
            obj.DEG2RAD = (pi / 180);
            obj.wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(obj.ecefFrame);
            obj.moon = CelestialBodyFactory.getMoon();
            obj.sun = CelestialBodyFactory.getSun();
            obj.groundStationsList = [];
            obj.numberOfGroundStation = 0;
            obj.spacecraftsList = [];
            obj.propagatorsList = [];
            tV = inputStructure.startTimeVector;
            obj.startTime = AbsoluteDate(tV(1), tV(2), tV(3), tV(4), tV(5),...
                tV(6), TimeScalesFactory.getUTC());
            obj.duration = inputStructure.OdDurationInMinutes / 60; % duration in hours
            obj.durationAfterEstimation = ...
                inputStructure.durationAfterEstimation / 60.0; % after estimation in hours
            obj.endEstimation = obj.startTime.shiftedBy(obj.duration * 3600.0); % final estimation time
            obj.endTime = obj.endEstimation.shiftedBy(obj.durationAfterEstimation * 3600.0); % final getting results time
            obj.plotStepTime = inputStructure.plotStepTime; % in sec
            obj.mass = inputStructure.satelliteMass;
            obj.referenceOrbit = KeplerianOrbit(...
                Constants.WGS84_EARTH_EQUATORIAL_RADIUS + ...
                (inputStructure.orbitAltitude * 1000), ...
                inputStructure.orbitEccentricity,  ...
                inputStructure.orbitInclination * obj.DEG2RAD,  ...
                inputStructure.orbitPregeeArgument * obj.DEG2RAD,  ...
                inputStructure.orbitRaan * obj.DEG2RAD,  ...
                inputStructure.orbitTrueAnomaly * obj.DEG2RAD,  ...
                PositionAngleType.TRUE, ...
                obj.eciFrame, ...
                obj.startTime,  ...
                Constants.WGS84_EARTH_MU);
            obj.nominalAltitudes = [inputStructure.orbitAltitude - 100,...
                inputStructure.orbitAltitude + 100];
            obj.noiseSeed = inputStructure.noiseSeed;
            obj.elevationMask = inputStructure.maskAngle * obj.DEG2RAD; % rad

             % Orbit propagator parameters
            obj.propagationMinTimeStep = inputStructure.propagationMinTimeStep;  % s
            obj.propagationMaxTimeStep = inputStructure.propagationMaxTimeStep;  % stimestep
            obj.vectorAbsoluteTolerance = inputStructure.vectorAbsoluteTolerance; % m
            obj.vectorRelativeTolerance = inputStructure.vectorRelativeTolerance; % m

            % other filter parameters
            obj.BLSestimatorConvergenceThreshold = inputStructure.BLSestimatorConvergenceThreshold;
            obj.BLSestimatorMaxIterations = inputStructure.BLSestimatorMaxIterations;
            obj.BLSestimatorMaxEvaluations = inputStructure.BLSestimatorMaxEvaluations;
            obj.showActiveSatellites = inputStructure.showActiveSatellites;

            obj.rangeMaxError = inputStructure.stationRangeMaxError;
            obj.rangeSigma = inputStructure.stationRangeSigma;
            obj.rangeWeight = inputStructure.stationRangeWeight;
            obj.stationSigma = inputStructure.stationSigma; % RangeRate sigma
            obj.stationBaseWeight = inputStructure.stationBaseWeight;
            obj.stationisTwoWay = inputStructure.stationisTwoWay;
            obj.stationRefractio = inputStructure.stationRefraction;
            obj.stationSampleTime = inputStructure.stationSampleTime;
            obj.stationStd = inputStructure.stationStd;

            obj.gpsSampleTime = inputStructure.gpsSampleTime;    % Observation parameter
            obj.gpsSigma = inputStructure.gpsSigma;         % filter parameter
            obj.gpsBaseWeight = inputStructure.gpsBaseWeight;    % filter parameter
            obj.gpsPositionStd = inputStructure.gpsPositionStd; % GPS Position Error ~ 3-sigma: 100 [m]
            obj.gpsVelocitystd = inputStructure.gpsVelocitystd; % GPS Velocity Error
            obj.GPSMaxBurstSize = 1;    % maximum number of selected dates in a burst
            obj.GPSHighRateStep = 1.0 * obj.gpsSampleTime; % step between two consecutive dates within a burst (s)
            obj.GPSBurstPeriod  = 1.0 * obj.gpsSampleTime; % period between the start of each burst (s)

            obj.isISL = inputStructure.isISL;
            obj.islSigma = inputStructure.islSigma;      % filter parameter
            obj.islBaseWeight = inputStructure.islBaseWeight;    % filter parameter
            obj.isTwoWay = inputStructure.isTwoWay;        % Observation parameter
            obj.islSampleTime = inputStructure.islSampleTime;   % Observation parameter
            obj.islScheduler = inputStructure.islScheduler; % alternative 'InView'
            obj.islStd = inputStructure.islStd;        % Environment parameter

            if strcmpi(inputStructure.islLink, 'G01')
                obj.islWaveLength = Frequency.G01.getWavelength();
            elseif strcmpi(inputStructure.islLink, 'G05')
                obj.islWaveLength = Frequency.G05.getWavelength();
            end


        end

        function satEstimator = getPodEstimator(obj, propagatorBuilder)
            import org.hipparchus.linear.*;
            import org.hipparchus.optim.nonlinear.vector.leastsquares.*;
            import org.orekit.estimation.leastsquares.*;
            import org.orekit.estimation.*;
            matrixDecomposer = QRDecomposer(1e-11);
            optimizer = GaussNewtonOptimizer(matrixDecomposer, false);
            satEstimator = BatchLSEstimator(optimizer, propagatorBuilder);
            satEstimator.setParametersConvergenceThreshold(obj.BLSestimatorConvergenceThreshold)
            satEstimator.setMaxIterations(obj.BLSestimatorMaxIterations)
            satEstimator.setMaxEvaluations(obj.BLSestimatorMaxEvaluations)
        end

        function feedPodMeasurements(obj, estimator, allOrbits, startTime, duration, ...
                measurementsType, linkList)
            if strcmp(measurementsType,"GPS")
                meas = obj.generateMeasurements(allOrbits, [], [], 0.0, ...
                    obj.elevationMask, "GPS",...
                    startTime, duration, obj.eciFrame, linkList);
                i = 0;
                for m = 1 : size(meas)
                    estimator.addMeasurement(meas(m));
                    i = i + 1;
                end
                fprintf("%1.0f measurements added to estimator \n",i);
            else
                for j = 1:obj.numberOfGroundStation
                    station = obj.groundStationsList(j);
                    stationName = obj.groundStationsList(j).getBaseFrame().getName();
                    altitude = obj.groundStationsList(j).getBaseFrame().getPoint().getAltitude();
                    meas = obj.generateMeasurements(allOrbits, station, stationName, altitude, ...
                        obj.elevationMask, measurementsType, ...
                        startTime, duration, obj.eciFrame, []);
                    i = 0;
                    for m = 1 : size(meas)
                        estimator.addMeasurement(meas(m));
                        i = i + 1;
                    end
                    fprintf("%1.0f measurements added to estimator \n",i);
                end

            end
        end

        function [orbits, selectedOrbits, selectedOrbitsIndexes, nominalLinksList] = ...
                buildISLStar(obj, num_plane, num_sat, F, activeSats)
            %** reference satellite must be assigned before calling this method!!!
            %inputs: number of plane: int
            %        number of satellites per plane: int
            %        walker phasing parameter: int (between zero and np-1)
            %        reference satellite: orekit KeplerianOrbit
            %        list of satellites with ISL link: list[dtype=int] length(between 1 and np*ns)
            if (num_plane-floor(num_plane)~=0)
                num_plane = floor(num_plane);
            end
            if (num_sat-floor(num_sat)~=0)
                num_sat = floor(num_sat);
            end
            if (F-floor(F)~=0)
                F = floor(F);
            end
            if F < 0
                F = 0;
            end
            if F > num_plane
                F = num_plane;
            end
            for i = 1:size(activeSats)
                if (activeSats(i)-floor(activeSats(i))~=0 )
                    activeSats(i)=floor(activeSats(i));
                end
                if  activeSats(i) < 0
                    activeSats(i) = 0;
                end
                if  activeSats(i) > num_plane*num_sat
                    activeSats(i) = num_plane*num_sat;
                end
            end
            activeSats = unique((activeSats));

            nominalLinksList = zeros(size(activeSats,2),5);
            allContributedSatellitesList = [];
            selectedOrbitsIndexes = [];
            selectedOrbits = [];
            orbits = [];
            refSat = obj.referenceOrbit;
            raan0 = refSat.getRightAscensionOfAscendingNode() * (180.0 / pi);
            ta0 = refSat.getTrueAnomaly() * (180.0 / pi);
            import org.orekit.utils.*;
            import org.orekit.orbits.*;
            for i = 1:num_plane
                for j = 1:num_sat
                    raan = raan0 + i * 180.0 / num_plane;
                    ta = ta0 + j * 360.0 / num_sat + i * 360 * F / (num_sat * num_plane);
                    ta = mod(ta,360.0);
                    if raan > 180.0
                        raan = raan - 180.0;
                    end
                    newOrbit = KeplerianOrbit(refSat.getA(), refSat.getE(), refSat.getI(), ...
                        refSat.getPerigeeArgument(), raan * (pi / 180.0), ...
                        ta * (pi / 180), PositionAngleType.TRUE, ...
                        refSat.getFrame(), refSat.getPVCoordinates().getDate(), ...
                        Constants.WGS84_EARTH_MU);
                    orbits = [orbits;newOrbit];
                end
            end
            fprintf(" ");
            fprintf("walker-star constellation created (%1.0f, %1.0f, %1.0f).\n",num_sat*num_plane, num_plane, F);

            for i = 1: size(activeSats,2)
                currentSatellite = activeSats(i);
                nominalLinksList(i , 1) = currentSatellite;
                % finding neigbhor satellites
                currentPlane = floor((currentSatellite - 1)/(num_sat) ) + 1;

                if (currentPlane == 1)
                    nominalLinksList(i, 2) = -1;
                    nominalLinksList(i, 4) = currentSatellite + num_sat;
                elseif (currentPlane == num_plane)
                    nominalLinksList(i, 2) = currentSatellite - num_sat;
                    nominalLinksList(i, 4) = -1;
                else
                    nominalLinksList(i, 2) = currentSatellite - num_sat;
                    nominalLinksList(i, 4) = currentSatellite + num_sat;
                end
                if  (currentSatellite == (currentPlane - 1)*num_sat + 1)
                    nominalLinksList(i, 3) = currentSatellite + 1;
                    nominalLinksList(i, 5) = currentSatellite + num_sat - 1;
                elseif (currentSatellite == currentPlane*num_sat)
                    nominalLinksList(i, 3) = currentSatellite - num_sat + 1;
                    nominalLinksList(i, 5) = currentSatellite - 1;
                else
                    nominalLinksList(i, 3) = currentSatellite + 1;
                    nominalLinksList(i, 5) = currentSatellite - 1;
                end
                fprintf(" ")
                fprintf("Satellite %1.0f is equipped with ISL links: (%1.0f,%1.0f,%1.0f,%1.0f)\n",currentSatellite , ...
                    nominalLinksList(i , 2), nominalLinksList(i , 3), nominalLinksList(i , 4), nominalLinksList(i , 5) )
                fprintf("         %1.0f         \n",nominalLinksList(i, 3))
                fprintf("          |          \n")
                fprintf("          |          \n")
                fprintf("%1.0f-------%1.0f---------%1.0f\n",nominalLinksList(i, 2), nominalLinksList(i, 1), nominalLinksList(i, 4))
                fprintf("          |          \n")
                fprintf("          |          \n")
                fprintf("         %1.0f          \n",nominalLinksList(i, 5))
            end
            for i = 1:size(activeSats,2)
                allContributedSatellitesList = [allContributedSatellitesList;activeSats(i)]; %#ok<AGROW>
                for j = 1:4
                    if nominalLinksList(i , j ) == -1
                        ...
                    else
                        allContributedSatellitesList = [allContributedSatellitesList; nominalLinksList(i, j + 1)]; %#ok<AGROW>
                    end
                end
            end
            uniqueIndexes = unique(allContributedSatellitesList);
            for i = 1 : size(uniqueIndexes,1)
                if  (uniqueIndexes(i) ~= -1)
                    selectedOrbits = [selectedOrbits;orbits(uniqueIndexes(i))];
                    selectedOrbitsIndexes = [selectedOrbitsIndexes;uniqueIndexes(i)];
                end
            end


        end

        function satPropagator = getPropagator(obj, propagatorType, initialOrbit, level)
            import org.orekit.time.*;
            import org.orekit.attitudes.*;
            import org.orekit.orbits.*;
            import org.orekit.propagation.*;
            import org.orekit.propagation.conversion.*;
            import org.hipparchus.ode.nonstiff.*;
            import org.orekit.propagation.numerical.*;
            import org.orekit.forces.gravity.potential.*;
            import org.orekit.forces.gravity.*;
            import org.orekit.models.earth.atmosphere.data.*;
            import org.orekit.models.earth.atmosphere.*;
            import org.orekit.forces.drag.*;
            import org.orekit.forces.radiation.*;
            import org.orekit.utils.*;
            minStep = obj.propagationMinTimeStep;
            maxStep = obj.propagationMaxTimeStep;
            vecAbsoluteTolerance = obj.vectorAbsoluteTolerance;
            vecRelativeTolerance = obj.vectorRelativeTolerance;
            nadirPointing = NadirPointing(obj.eciFrame, obj.wgs84Ellipsoid);
            initialCartesianOrbit = CartesianOrbit(SpacecraftState(initialOrbit, obj.mass).getPVCoordinates(obj.eciFrame),...
                obj.eciFrame, obj.wgs84Ellipsoid.getGM());

            if strcmpi(propagatorType,'Builder')
                integratorBuilder = DormandPrince853IntegratorBuilder(minStep, maxStep, vecAbsoluteTolerance);
                satPropagator = NumericalPropagatorBuilder(initialCartesianOrbit, integratorBuilder, PositionAngleType.TRUE,...
                    1.0);
                satPropagator.setAttitudeProvider(nadirPointing);
                satPropagator.setMass(obj.mass);
            else
                thisintegrator = DormandPrince853Integrator(minStep, ...
                    maxStep, ...
                    vecAbsoluteTolerance, ...
                    vecRelativeTolerance);
                satPropagator = NumericalPropagator(thisintegrator);
                satPropagator.setInitialState(SpacecraftState(initialCartesianOrbit, obj.mass));
                satPropagator.setAttitudeProvider(nadirPointing);
            end
            % determine the level of the propagator
            if strcmpi(level,"Low") || strcmpi(level,"l")
                propagatorCase = 1;
            elseif strcmpi(level,"Medium") || strcmpi(level,"m")
                propagatorCase = 2;
            elseif strcmpi(level,"High") || strcmpi(level,"H")
                propagatorCase = 3;
            else
                propagatorCase = 0;
            end
            if (propagatorCase == 0)
                gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(1, 1);
                gravityAttractionModel = HolmesFeatherstoneAttractionModel(obj.ecefFrame, gravityProvider);
                satPropagator.addForceModel(gravityAttractionModel);
            elseif (propagatorCase == 1)
                gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(2, 2, obj.startTime);
                gravityAttractionModel = HolmesFeatherstoneAttractionModel(obj.ecefFrame, gravityProvider);
                satPropagator.addForceModel(gravityAttractionModel);
            elseif (propagatorCase == 2)
                gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(21, 21, obj.startTime);
                gravityAttractionModel = HolmesFeatherstoneAttractionModel(obj.ecefFrame, gravityProvider);
                % Atmospheric drag
                AVERAGE = javaMethod('valueOf','org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation$StrengthLevel','AVERAGE');
                msafe = MarshallSolarActivityFutureEstimation( MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES,AVERAGE);
                atmosphere = DTM2000(msafe, obj.sun, obj.wgs84Ellipsoid);
                isotropicDrag = IsotropicDrag(0.02, 2.2);
                dragForce = DragForce(atmosphere, isotropicDrag);
                % Solar radiation pressure
                isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(0.02, 1.0);
                solarRadiationPressure = SolarRadiationPressure(obj.sun, obj.wgs84Ellipsoid, ...
                    isotropicRadiationSingleCoeff);
                % Third body attraction model
                moon_3dbodyattraction = ThirdBodyAttraction(obj.moon);
                sun_3dbodyattraction = ThirdBodyAttraction(obj.sun);
                satPropagator.addForceModel(gravityAttractionModel)
                satPropagator.addForceModel(dragForce)
                satPropagator.addForceModel(solarRadiationPressure)
                satPropagator.addForceModel(moon_3dbodyattraction)
                satPropagator.addForceModel(sun_3dbodyattraction)

            elseif (propagatorCase == 3)
                % Earth gravity field with degree 64 and order 64
                gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(64, 64, obj.startTime);
                gravityAttractionModel = HolmesFeatherstoneAttractionModel(obj.ecefFrame, gravityProvider);
                % Third body attraction model
                moon_3dbodyattraction = ThirdBodyAttraction(obj.moon);
                sun_3dbodyattraction = ThirdBodyAttraction(obj.sun);

                % Solar radiation pressure
                isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(0.02, 1.0);
                solarRadiationPressure = SolarRadiationPressure(obj.sun, obj.wgs84Ellipsoid, isotropicRadiationSingleCoeff);

                % Relativity
                relativity = Relativity(Constants.EIGEN5C_EARTH_MU);

                oceanicTides = OceanTides(obj.ecefFrame, ...
                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_MU, ...
                    5, 5, IERSConventions.IERS_2010, ...
                    TimeScalesFactory.getUT1(IERSConventions.IERS_2010, true));
                solidTidess = SolidTides(obj.ecefFrame, ...
                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_MU, ...
                    gravityProvider.getTideSystem(), ...
                    IERSConventions.IERS_2010, ...
                    TimeScalesFactory.getUT1(IERSConventions.IERS_2010, true), ...
                    [obj.sun, obj.moon]);

                % Atmospheric drag
                % from org.orekit.models.earth.atmosphere import NRLMSISE00
                % atmosphere = NRLMSISE00(msafe, sun, wgs84Ellipsoid)
                AVERAGE = javaMethod('valueOf','org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation$StrengthLevel','AVERAGE');
                msafe = MarshallSolarActivityFutureEstimation( ...
                    MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES, AVERAGE);
                atmosphere = DTM2000(msafe, obj.sun, obj.wgs84Ellipsoid);
                isotropicDrag = IsotropicDrag(0.02, 2.2);
                dragForce = DragForce(atmosphere, isotropicDrag);
                satPropagator.addForceModel(gravityAttractionModel);
                satPropagator.addForceModel(moon_3dbodyattraction);
                satPropagator.addForceModel(sun_3dbodyattraction);
                satPropagator.addForceModel(solarRadiationPressure);
                satPropagator.addForceModel(relativity);
                satPropagator.addForceModel(dragForce);
                satPropagator.addForceModel(oceanicTides);
                satPropagator.addForceModel(solidTidess);
            end

        end

        function feedIslMeasurements(obj, estimator, underODPropagator, linkedpropagators, ...
                startTime, duration, islType, linkList)
            import org.orekit.time.*;
            import org.hipparchus.random.*;
            import org.hipparchus.linear.*;
            import org.orekit.estimation.measurements.generation.*;
            import org.orekit.gnss.*;
            import org.orekit.propagation.events.*;
            initIslTime = startTime;
            finalIslTime = initIslTime.shiftedBy(duration * 3600.0);
            small = 1e-18;
            seed = obj.noiseSeed;
            random_generator = Well19937a(seed);
            gaussian_generator = GaussianRandomGenerator(random_generator);
            covariance = MatrixUtils.createRealDiagonalMatrix([obj.islStd^2, obj.islStd^2]);
            noise_source = CorrelatedRandomVectorGenerator(covariance, small, gaussian_generator);
            gatherer = GatheringSubscriber();
            generator = Generator();
            generator.addSubscriber(gatherer);
            local = underODPropagator;
            localEphemerisGenerator = local.getEphemerisGenerator();
            local.propagate(initIslTime , finalIslTime);
            localBoundedPropagator = localEphemerisGenerator.getGeneratedEphemeris();
            satellite = generator.addPropagator(localBoundedPropagator);
            Rsatellites = [];
            for i = 1 : size(linkedpropagators, 1)
                remote = linkedpropagators(i);
                remoteEphemerisGenerator = remote.getEphemerisGenerator();
                remote.propagate(initIslTime , finalIslTime);
                remoteBoundedPropagator = remoteEphemerisGenerator.getGeneratedEphemeris();
                Rsatellites = [Rsatellites; generator.addPropagator(remoteBoundedPropagator)]; %#ok<AGROW>
            end
            for j = 1 : size(Rsatellites, 1)
                fixed_step_selector = FixedStepSelector(obj.islSampleTime , TimeScalesFactory.getUTC());
                if strcmpi(islType, 'Range')
                    Builder = InterSatellitesRangeBuilder(noise_source, satellite, ...
                        Rsatellites(j,1), obj.isTwoWay, ...
                        obj.islSigma, obj.islBaseWeight);
                else
                    Builder = InterSatellitesPhaseBuilder(noise_source, satellite, ...
                        Rsatellites(j), obj.islWaveLength, ...
                        obj.islSigma, obj.islBaseWeight);
                end
                if strcmpi(obj.islScheduler,'InView')
                    interDetector = InterSatDirectViewDetector( obj.wgs84Ellipsoid , linkedpropagators(j));
                    scheduler = EventBasedScheduler(Builder, fixed_step_selector, underODPropagator, ...
                        interDetector, ...
                        SignSemantic.FEASIBLE_MEASUREMENT_WHEN_POSITIVE);
                    generator.addScheduler(scheduler);
                else
                    generator.addScheduler(ContinuousScheduler(Builder, fixed_step_selector));
                end
                generator.generate(initIslTime, finalIslTime);
                measurements = gatherer.getGeneratedMeasurements();
                measurementIterator = measurements.iterator();
                for mn = 1:measurements.size()
                    measObject = measurementIterator.next();
                    estimator.addMeasurement(measObject)

                end
                fprintf("%1.0f measurements feeded to satellite #%1.0f estimator from ISL %s observations of satellite #%1.0f \n",mn,...
                    linkList(1), islType, linkList(j+1))
                clear Builder;
            end
        end

        function measurementslist =  generateMeasurements(obj, allOrbits, station, station_name, ...
                altitude, elevAngle, meas_type, t0, duration, eciFrame, linkList)
            import org.orekit.time.*;
            import org.hipparchus.random.*;
            import org.hipparchus.linear.*;
            import org.orekit.estimation.measurements.generation.*;
            import org.orekit.gnss.*;
            import org.orekit.propagation.events.*;
            import org.orekit.models.earth.*;
            import org.orekit.propagation.events.handlers.*;
            two_way = obj.stationisTwoWay;
            withRefraction = obj.stationRefractio;
            step = obj.stationSampleTime;
            seed = obj.noiseSeed;
            tf = t0.shiftedBy(3600.0 * duration);
            measurementslist = [];
            for proIndex = 1:size(allOrbits,1)
                gatherer = GatheringSubscriber();
                generator = Generator();
                generator.addSubscriber(gatherer);
                for proIndex2 = 1:size(allOrbits,1)
                    if proIndex2 == proIndex
                        theODp = obj.getPropagator("Get", allOrbits(proIndex2), "H");
                        satellite = generator.addPropagator(theODp);
                    else
                        generator.addPropagator(obj.getPropagator("Get", allOrbits(proIndex2), "H"));
                    end
                end
                fixed_step_selector = FixedStepSelector(step, TimeScalesFactory.getUTC());
                sigma = obj.stationSigma;
                base_weight = obj.stationBaseWeight;
                small = 1e-10;
                random_generator = Well19937a(seed);
                gaussian_generator = GaussianRandomGenerator(random_generator);
                covariance = MatrixUtils.createRealDiagonalMatrix(...
                    [(obj.stationStd * obj.stationStd), ...
                    (obj.stationStd * obj.stationStd)]);
                noise_source = CorrelatedRandomVectorGenerator(covariance, small, gaussian_generator);
                if strcmpi(meas_type,'RANGE')
                    builder = RangeBuilder(noise_source, station, two_way, sigma, base_weight, satellite);
                elseif strcmpi(meas_type,'RANGERATE')
                    builder = RangeRateBuilder(noise_source, station, two_way, sigma, base_weight, satellite);
                elseif strcmpi(meas_type,'AZEL')
                    sigmaAzEl = JArray_double([sigma, sigma]);
                    weightAzEl = JArray_double([base_weight, base_weight]);
                    builder = AngularAzElBuilder(noise_source, station, sigmaAzEl, weightAzEl, satellite);
                elseif strcmpi(meas_type,'RADEC')
                    sigmaRaDec = JArray_double([sigma, sigma]);
                    weightRaDec = JArray_double([base_weight, base_weight]);
                    builder = AngularRaDecBuilder(noise_source, station, eciFrame, sigmaRaDec, weightRaDec, satellite);
                elseif strcmpi(meas_type,'GPS')
                    sigmaP = obj.gpsPositionStd ^ 2;
                    sigmaV = obj.gpsVelocitystd ^ 2;
                    GPSCov = MatrixUtils.createRealDiagonalMatrix([sigmaP, sigmaP, sigmaP, sigmaV, sigmaV, sigmaV]);
                    GPSNoise = CorrelatedRandomVectorGenerator(GPSCov, (small), gaussian_generator);
                    GPSBuilder = PVBuilder(GPSNoise, obj.gpsSigma, obj.gpsSigma*0.01, obj.gpsBaseWeight, satellite);
                    GPSsample = BurstSelector(obj.GPSMaxBurstSize, obj.GPSHighRateStep, obj.GPSBurstPeriod, TimeScalesFactory.getUTC());
                    GPSScheduler = ContinuousScheduler(GPSBuilder, GPSsample);
                end
                if strcmpi(meas_type,'GPS')
                    generator.addScheduler(GPSScheduler)
                else
                    if withRefraction
                        elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(...
                            elevAngle).withRefraction(EarthITU453AtmosphereRefraction(altitude)).withHandler(ContinueOnEvent());
                    else
                        elevation_detector = ElevationDetector(station.getBaseFrame()).withConstantElevation(elevAngle).withHandler(...
                            ContinueOnEvent());
                    end
                    scheduler = EventBasedScheduler(builder, fixed_step_selector, ...
                        theODp, ...
                        elevation_detector,...
                        SignSemantic.FEASIBLE_MEASUREMENT_WHEN_POSITIVE);
                    generator.addScheduler(scheduler);
                end
                generator.generate(t0, tf);
                measurements = gatherer.getGeneratedMeasurements();
                measurementIterator = measurements.iterator();
                measTime = [];
                for j = 1 : measurements.size()
                    meas = measurementIterator.next();
                    measurementslist = [measurementslist; meas];
                    if j == 1
                        measTime = [measTime; meas];
                    end
                end
                try
                    if strcmpi(meas_type,'GPS')
                        measTime = [measTime; meas];
                        fprintf("%1.0f number of GPS observations are generated for satellite #%1.0f, first at: %s and last at: %s \n",j, ...
                            linkList(proIndex), measTime(1).getDate().toString(), measTime(2).getDate().toString());
                    else
                        measTime = [measTime; meas];
                        fprintf("%1.0f number of %s observations are generated for satellite #%1.0f, first at: %s and last at: %s, - station %s \n",j, ...
                            meas_type, linkList(proIndex), measTime(1).getDate().toString(), measTime(2).getDate().toString(), station_name);
                    end
                catch
                    ...
                end
            end
        end

        function positionResidual = getRICPositionDifference(obj, prop1, prop2)
            import org.orekit.frames.*;
            wholeDuration = obj.endTime.durationFrom(obj.startTime);
            currentDateTime = obj.startTime;
            positionResidual = zeros(floor(wholeDuration/obj.plotStepTime) + 1,4);
            mn = 0;
            while currentDateTime.compareTo(obj.endTime) <= 0
                mn = mn + 1;
                truePosition = prop1.propagate(currentDateTime).getPVCoordinates().getPosition();
                estimatedPosition = prop2.propagate(currentDateTime).getPVCoordinates().getPosition();
                positionDifference = truePosition.subtract(estimatedPosition);
                errorInLVLH = LOFType.LVLH.rotationFromInertial(prop1.propagate(currentDateTime).getPVCoordinates()).applyTo(positionDifference);
                errorVector = [errorInLVLH.getX(), errorInLVLH.getY(), errorInLVLH.getZ(), errorInLVLH.getNorm()];
                positionResidual(mn,:) =  errorVector;
                currentDateTime = currentDateTime.shiftedBy(obj.plotStepTime);
            end
        end

        function additionalFigureInfo = getFigureTitle(obj, linkNumber)
            % if obj.isISL
            %     additionalFigureInfo = ['ISL is ON, ISL sigma = '  num2str(obj.islSigma)  ...
            %         ', ISL BaseWeight = '  num2str(obj.islBaseWeight)  ', GPS sigma = '  num2str(obj.gpsSigma)  ...
            %         ', GPS BaseWeight = '  num2str(obj.gpsBaseWeight)  ', GPS noise sigma = '  num2str(obj.gpsPositionStd)  ...
            %         ',  for Sat:#'   num2str(linkNumber)];
            % else
            %     additionalFigureInfo = ['ISL is OFF, GPS sigma = '  num2str(obj.gpsSigma) ...
            %         ', GPS BaseWeight = '  num2str(obj.gpsBaseWeight)  ' GPS noise sigma = '  num2str(obj.gpsPositionStd) ...
            %         ',  for Sat:#'   num2str(linkNumber)];
            % end

            if obj.isISL
                additionalFigureInfo = ['ISL is ON,   Error History for Sat:#'   num2str(linkNumber) ' in RIC Frame'];
            else
                additionalFigureInfo = ['ISL is OFF,  Error History for Sat:#'   num2str(linkNumber) ' in RIC Frame'];
            end

        end

        function displayActiveSatellitesConfiguration(obj, activeSatellites, selectedOrbitsIndexes, linkList)
            %earhtMarbleFigure = self.getBlueMarbleFigure('earth.jpeg')
            epochECI2ECEF = obj.eciFrame.getTransformTo(obj.ecefFrame, obj.startTime);
            mn = 0;
            for i  = 1:size(linkList, 2)
                if selectedOrbitsIndexes(i) == linkList(1)
                    mainSat = i;
                end
                if selectedOrbitsIndexes(i) == -1
                    mn = mn +1;
                end
            end
            mainSatEci = activeSatellites(mainSat).getPVCoordinates().getPosition();
            mainSatECEF = epochECI2ECEF.getRotation().applyTo(mainSatEci);
            allPositionECEF = mainSatECEF;
            for j = 1:size(activeSatellites,1)
                if j == mainSat
                    ...
                else
                    eciPos = activeSatellites(j).getPVCoordinates().getPosition();
                    truePositionECEF = epochECI2ECEF.getRotation().applyTo(eciPos);
                    allPositionECEF = [allPositionECEF; truePositionECEF];
                end
            end
            geoMap = imread("blue.PNG");
            figure()
            [X,Y,Z] = sphere(50);
            hold on
            rotate3d
            axis equal
            mesh(X.*6378137,Y.*6378137,Z.*6378137,flipud(geoMap),'FaceColor','texturemap',EdgeColor='none')
            plot3(allPositionECEF(1).getX(), allPositionECEF(1).getY(), allPositionECEF(1).getZ(), 'r.','MarkerSize',25)
            disss = 400000;
            text(allPositionECEF(1).getX()+ disss, allPositionECEF(1).getY()+ disss, allPositionECEF(1).getZ() + disss, ...
                num2str(linkList(1)),'Color','r', 'FontSize',12)
            for i = 2:size(allPositionECEF)
                plot3(allPositionECEF(i).getX(), allPositionECEF(i).getY(), allPositionECEF(i).getZ(), 'k.','MarkerSize',25)
                text(allPositionECEF(i).getX()+ disss, allPositionECEF(i).getY()+ disss,...
                    allPositionECEF(i).getZ() + disss,num2str(linkList(i)), 'FontSize',12)
            end
            axis off
        end

        function plotReferenceOrbitTrajectory(obj, durationInHours, plotSampleTime)

            currentDateTime = obj.startTime;
            thisEnd = currentDateTime.shiftedBy(durationInHours*3600.0);
            wholeDuration = thisEnd.durationFrom(currentDateTime);
            positionArrayECI = zeros(floor(wholeDuration/plotSampleTime) + 1,3);
            positionArrayECEF = zeros(floor(wholeDuration/plotSampleTime) + 1,3);
            thisOrbit = obj.referenceOrbit;
            thisPropagator = obj.getPropagator('Get', thisOrbit, 'H');
            thisEphemerisGenerator = thisPropagator.getEphemerisGenerator();
            thisPropagator.propagate(currentDateTime, thisEnd);
            thisBoundedPropagator = thisEphemerisGenerator.getGeneratedEphemeris();
            mn = 0;
            while currentDateTime.compareTo(thisEnd) <= 0
                mn = mn + 1;
                thisECI2ECEF = obj.eciFrame.getTransformTo(obj.ecefFrame, currentDateTime);
                thisPositionECI = thisBoundedPropagator.propagate(currentDateTime).getPVCoordinates().getPosition();
                thisPositionECEF = thisECI2ECEF.getRotation().applyTo(thisPositionECI);
                currentDateTime = currentDateTime.shiftedBy(plotSampleTime);
                positionArrayECI(mn,1) = thisPositionECI.getX();
                positionArrayECI(mn,2) = thisPositionECI.getY();
                positionArrayECI(mn,3) = thisPositionECI.getZ();
                positionArrayECEF(mn,1) = thisPositionECEF.getX();
                positionArrayECEF(mn,2) = thisPositionECEF.getY();
                positionArrayECEF(mn,3) = thisPositionECEF.getZ();
            end
            geoMap = imread("blue.PNG");
            figure()
            [X,Y,Z] = sphere(50);
            hold on
            rotate3d
            axis equal
            mesh(X.*6378137,Y.*6378137,Z.*6378137,flipud(geoMap),'FaceColor','texturemap',EdgeColor='none')
            plot3(positionArrayECEF(:,1),positionArrayECEF(:,2),positionArrayECEF(:,3), 'LineWidth',3, 'Color','k')
            axis off
        end

        function displayConstellationConfiguration(obj, orbits)
            %earhtMarbleFigure = self.getBlueMarbleFigure('earth.jpeg')
            epochECI2ECEF = obj.eciFrame.getTransformTo(obj.ecefFrame, obj.startTime);
            allPositionECEF = zeros(size(orbits,1),3);
            for j = 1:size(orbits,1)
                eciPos = orbits(j).getPVCoordinates().getPosition();
                thisPositionECEF = epochECI2ECEF.getRotation().applyTo(eciPos);
                allPositionECEF(j,1) = thisPositionECEF.getX();
                allPositionECEF(j,2) = thisPositionECEF.getY();
                allPositionECEF(j,3) = thisPositionECEF.getZ();
            end
            geoMap = imread("blue.PNG");
            figure()
            [X,Y,Z] = sphere(50);
            hold on
            rotate3d
            axis equal
            mesh(X.*6378137,Y.*6378137,Z.*6378137,flipud(geoMap),'FaceColor','texturemap',EdgeColor='none')
            for i = 1:size(allPositionECEF,1)
                plot3(allPositionECEF(i,1), allPositionECEF(i,2), allPositionECEF(i,3), 'k.','MarkerSize',16)
            end
            axis off
        end

        function thisFigure = getErrorFigure(obj, ricposdiff, figureInfo)
            allSteps = (1:size(ricposdiff(:,4))).*(obj.plotStepTime / 60.0);
            thisFigure = figure();
            subplot(1,4,1:3)
            plot(allSteps, ricposdiff(:,4), 'k', 'LineWidth',5)
            hold on
            plot(allSteps, ricposdiff(:,1),'LineStyle','--', 'LineWidth',3)
            plot(allSteps, ricposdiff(:,2),'LineStyle','--', 'LineWidth',3)
            plot(allSteps, ricposdiff(:,3),'LineStyle','--', 'LineWidth',3)
            grid on
            grid minor
            title(figureInfo)
            ylabel('Error in [m]')
            xlabel('Time lapsed from the epoch [min]')
            patch([0, 0, obj.duration * 60.0, obj.duration * 60.0],...
                [gca().YLim(1), gca().YLim(2), gca().YLim(2), gca().YLim(1)],'g', 'edgecolor','none')
            alpha(0.4)
            legend({'Cumulative Error', 'In-tack Error', 'Cross-track Error', 'radial Error'})
            set(gca, 'FontName', 'Times New Roman')
            set(gca, 'FontSize', 20)
            subplot(1,4,4)
            histogram(ricposdiff(:,4))
            xlabel('Error value in [m]')
            ylabel('No. of Samples')
            title(['Norm Error RMS: ' num2str(rms(ricposdiff(:,4)))])
            set(gca, 'FontName', 'Times New Roman')
            set(gca, 'FontSize', 20)
        end

        function [meshGrid, EcefToEnuCells]= generateMeshGrid(obj, initialNumberOfPoints)
            % Licensed
            PHI = (1 + sqrt(5)) / 2;
            NUM_TRIANGLES_IN_AN_ICOSAHEDRON = 20;
            NUM_DIMENSIONS = 3;
            NUM_VERTICES_IN_A_TRIANGLE = 3;
            NUM_NEW_TRIANGLES_PER_TRIANGLE = 4;
            VERTICES = [0, PHI, 1; 0, -PHI, 1; 0, PHI, -1; 0, -PHI, -1; 1, 0, PHI; ...
                -1, 0, PHI; 1, 0, -PHI; -1, 0, -PHI; PHI, 1, 0; -PHI, 1, 0; PHI, -1, ...
                0; -PHI, -1, 0] / norm([1, PHI]);
            As = [2; 5; 9; 7; 11; 4; 6; 2; 1; 5; 3; 9; 8; 7; 12; 12; 6; 1; 3; 10];
            Bs = [4; 11; 5; 9; 7; 2; 12; 5; 6; 9; 1; 7; 3; 4; 8; 6; 1; 3; 10; 8];
            Cs = [11; 2; 11; 11; 4; 12; 2; 6; 5; 1; 9; 3; 7; 8; 4; 10; 10; 10; 8; 12];
            triangles = reshape([VERTICES(As, :), VERTICES(Bs, :), VERTICES(Cs, :)], ...
                NUM_TRIANGLES_IN_AN_ICOSAHEDRON, NUM_DIMENSIONS, ...
                NUM_VERTICES_IN_A_TRIANGLE);
            numDivisions = ceil(obj.thisLogarithm((initialNumberOfPoints - 2) / 25, ...
                NUM_NEW_TRIANGLES_PER_TRIANGLE));
            for i = 1:numDivisions
                triangles = obj.divideTriangles(triangles, NUM_DIMENSIONS, ...
                    NUM_VERTICES_IN_A_TRIANGLE, NUM_NEW_TRIANGLES_PER_TRIANGLE);
            end
            [latGridInDegrees, longGridInDegrees] = ...
                obj.cartesianToSpherical(triangles(:, 1, 1), ...
                triangles(:, 2, 1), triangles(:, 3, 1));
            [sortedOnes, sortIndices1] = sort(latGridInDegrees);
            sortedTwos = longGridInDegrees(sortIndices1);
            uniqueOneIndices = obj.gridDistinctFromPrevious(sortedOnes);
            numUniqueOnes = numel(uniqueOneIndices);
            for i = 1:(numUniqueOnes - 1)
                cluster = uniqueOneIndices(i):(uniqueOneIndices(i + 1) - 1);
                [sortedTwos(cluster), sortIndices] = sort(sortedTwos(cluster));
                sortIndicesAsInts = uint32(sortIndices);
                sortedOnes(cluster) = sortedOnes(uniqueOneIndices(i) + ...
                    sortIndicesAsInts - 1);
            end
            finalCluster = uniqueOneIndices(numUniqueOnes):numel(sortedTwos);
            [sortedTwos(finalCluster), sortIndices] = sort(sortedTwos(finalCluster));
            sortIndicesAsInts = uint32(sortIndices);
            sortedOnes(finalCluster) = sortedOnes(uniqueOneIndices(numUniqueOnes) ...
                + sortIndicesAsInts - 1);
            uniqueTwoIndices = obj.gridDistinctFromPrevious(sortedTwos);
            uniquePairIndicesWithDuplicates = ...
                sort([uniqueOneIndices; uniqueTwoIndices]);
            uniquePairIndices = uniquePairIndicesWithDuplicates(...
                obj.gridDistinctFromPrevious(uniquePairIndicesWithDuplicates));
            prunedOnes = sortedOnes(uniquePairIndices);
            prunedTwos = sortedTwos(uniquePairIndices);
            meshGrid = zeros(size(prunedOnes,1),3);
            EcefToEnuCells = cell( length(prunedOnes) , 1);
            for i = 1:size(prunedOnes,1)
                meshGrid(i,:)  = obj.llaToEcef([prunedOnes(i) prunedTwos(i) 0.0]);
                lat = prunedOnes(i)*pi/180;
                lon = prunedTwos(i)*pi/180;
                EcefToEnuCells{i} = [
                    -sin(lon)         ,  cos(lon)         , 0;
                    -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat);
                    cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)];
            end
        end

        function propagators = getConstellationPropagators(obj, constellationOnEpoch, wholeDuration)
            startDate = constellationOnEpoch(1).getDate();
            endDate = startDate.shiftedBy(wholeDuration);
            propagators = [];
            for i = 1:size(constellationOnEpoch, 1)
                thisPropagator = obj.getPropagator('Get', constellationOnEpoch(i), 'H');
                thisEphemerisGenerator = thisPropagator.getEphemerisGenerator();
                thisPropagator.propagate(startDate, endDate);
                propagators =  [propagators; thisEphemerisGenerator.getGeneratedEphemeris()];
            end
        end

        function constellationStatesAtNewEpoch = getConstellaionStatesAt(obj, propagators, durationInSeconds)
            constellationStatesAtNewEpoch = [];
            oldDateTime = propagators(1).getMinDate();
            newDateTime = oldDateTime.shiftedBy(durationInSeconds);
            for i = 1:size(propagators, 1)
                constellationStatesAtNewEpoch = ...
                    [constellationStatesAtNewEpoch; ...
                    propagators(i).propagate(newDateTime)];
            end
        end

        function allPositionECEF =  getConstellationEcef(obj, orbits)
            %earhtMarbleFigure = self.getBlueMarbleFigure('earth.jpeg')
            epochECI2ECEF = obj.eciFrame.getTransformTo(obj.ecefFrame, orbits(1).getDate());
            allPositionECEF = zeros(size(orbits,1),3);
            for j = 1:size(orbits,1)
                eciPos = orbits(j).getPVCoordinates().getPosition();
                thisPositionECEF = epochECI2ECEF.getRotation().applyTo(eciPos);
                allPositionECEF(j,1) = thisPositionECEF.getX();
                allPositionECEF(j,2) = thisPositionECEF.getY();
                allPositionECEF(j,3) = thisPositionECEF.getZ();
            end

        end

        function LinkedSatellites = getLinkedSatellitesEcef(obj, constellationEcef, pointLlaDegrees, maskInDegrees)
            LinkedSatellites = [];
            mn = 0;
            pointEcef = obj.llaToEcef(pointLlaDegrees);
            for i = 1:size(constellationEcef, 1)
                [~, El, ~, ~] = obj.getAzEl(pointLlaDegrees, constellationEcef(i,:));
                if (El * (180/pi)) > maskInDegrees
                    index = i;
                    mn = mn + 1;
                    LinkedSatellites(mn,1:4) = [index, constellationEcef(i,:)];
                    LinkedSatellites(mn,5) = norm(pointEcef - constellationEcef(i,:));
                end
            end
        end

        function gmst = utc2gmst(obj, date1)
            % Compute the Julian date for the given input date vector.
            JD = juliandate(date1);

            % Compute UT1.
            UT1 = (JD-2451545.0)/36525;

            % Compute the Greenwich Mean Sidereal Time (GMST) [seconds].
            gmst = 67310.54841 + (876600*3600 + 8640184.812866)*UT1 + 0.093104*UT1^2 - 6.2e-6*UT1^3;

            % Convert GMST to radians and put in the range [0 2*pi].
            gmst = mod((gmst/240)*pi/180,2*pi);
        end

        function [beta, el, rho, R_ENU] = getAzEl(obj, stationLla, R_sat)
            lat = stationLla(1) * (pi/180);
            long = stationLla(2) * (pi/180);
            h_ellp = stationLla(3);


            % compute the position vector of the site in the ECEF frame
            R_site  = obj.llaToEcef(stationLla);

            % compute the range and range rate vectors in the ECEF frame
            rho_ECEF  = R_sat - R_site;
            %drho_ECEF = V_sat;

            % define the tranformation to the south east and up (SEZ) coordinate system
            mat(1,1) =  sin(lat)*cos(long);
            mat(1,2) =  sin(lat)*sin(long);
            mat(1,3) = -cos(lat);

            mat(2,1) = -sin(long);
            mat(2,2) =  cos(long);
            mat(2,3) =  0;

            mat(3,1) =  cos(lat)*cos(long);
            mat(3,2) =  cos(lat)*sin(long);
            mat(3,3) =  sin(lat);

            % transform to the SEZ coordinate system
            rho_SEZ  = mat*rho_ECEF';
            %drho_SEZ = mat*drho_ECEF';
            R_SEZ = rho_SEZ;
            R_ENU = [R_SEZ(2), -R_SEZ(1), R_SEZ(3)];

            % compute the range and range rate
            rho  = norm(rho_SEZ);
            %drho = dot(rho_SEZ,drho_SEZ)/rho;

            % compute the elevation angle
            el  = asin(rho_SEZ(3)/rho);

            temp  = sqrt( rho_SEZ(1)^2 + rho_SEZ(2)^2 );
            %dtemp = sqrt( drho_SEZ(1)^2 + drho_SEZ(2)^2 );

            %del = (drho_SEZ(3) - drho *sin(el)) / temp;

            % compute the azimuth angle
            if el == pi/2

                sinb =  drho_SEZ(2)/dtemp;
                cosb = -drho_SEZ(1)/dtemp;
                beta =  atan2(sinb,cosb);


            else

                sinb =  rho_SEZ(2)/temp;
                cosb = -rho_SEZ(1)/temp;
                beta =  atan2(sinb,cosb);

            end

            if beta < 0
                beta = beta + 2*pi;
            end

            %dbeta = (drho_SEZ(1)*rho_SEZ(2)-drho_SEZ(2)*rho_SEZ(1))/temp^2;

        end

        function [GDOP, PDOP, TDOP, HDOP, VDOP] = getGlobeDops(obj, ...
                initialNoPoints, starConstellation, el_mask, isFigureOn)
            sats = obj.getConstellationEcef(starConstellation);
            [sites, enus]= obj.generateMeshGrid(initialNoPoints);
            NoSites = size(sites,1);
            GDOP = zeros(NoSites,1);
            PDOP= size(NoSites,1);
            TDOP= size(NoSites,1);
            HDOP= size(NoSites,1);
            VDOP= size(NoSites,1);
            for i = 1:NoSites
                [GDOPt, PDOPt, TDOPt, HDOPt, VDOPt] = obj.getSiteDops(sites(i,:), enus{i}, sats, el_mask);
                GDOP(i,1) = GDOPt; PDOP(i,1) = PDOPt;
                TDOP(i,1) = TDOPt; HDOP(i,1) = HDOPt;
                VDOP(i,1) = VDOPt;
            end
            if isFigureOn
                figure()
                scatter3(sites(:,1),sites(:,2),sites(:,3),ones(size(sites,1),1).*10, GDOP)
                axis equal
                axis off
                cb = colorbar();
                ylabel(cb,'GDOP Value','FontName','Times New Roman','FontSize',16,'Rotation',270)
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 16)
            end
        end

        function [GDOP, PDOP, TDOP, HDOP, VDOP] = getGlobeDopsHistory(obj,...
                initialNoPoints, starConstellation, el_mask, durationInHours, sampleTimeInSeconds)
            linkedPropagators = [];
            for i = 1:size(starConstellation,1)
                thisPropagator = obj.getPropagator("Get", starConstellation(i), "L");
                thisGenerator = thisPropagator.getEphemerisGenerator();
                thisPropagator.propagate(obj.startTime, obj.startTime.shiftedBy(durationInHours*3600.0));
                linkedPropagators = [linkedPropagators; thisGenerator.getGeneratedEphemeris()];
            end
            [sites, enus]= obj.generateMeshGrid(initialNoPoints);
            NoTimeSamples = floor(durationInHours*3600/sampleTimeInSeconds);
            NoSites = size(sites,1);
            GDOP = zeros(NoSites,NoTimeSamples);
            PDOP= size(NoSites,NoTimeSamples);
            TDOP= size(NoSites,NoTimeSamples);
            HDOP= size(NoSites,NoTimeSamples);
            VDOP= size(NoSites,NoTimeSamples);
            for j = 1:NoTimeSamples
                sats = [];
                for z = 1:size(linkedPropagators,1)
                    sats = [sats; linkedPropagators(z).propagate(obj.startTime.shiftedBy(j*sampleTimeInSeconds))];
                end
                ecefs = obj.getConstellationEcef(sats);
                for i = 1:NoSites
                    [GDOPt, PDOPt, TDOPt, HDOPt, VDOPt] = obj.getSiteDops(sites(i,:), enus{i}, ecefs, el_mask);
                    GDOP(i,j) = GDOPt; PDOP(i,j) = PDOPt;
                    TDOP(i,j) = TDOPt; HDOP(i,j) = HDOPt;
                    VDOP(i,j) = VDOPt;
                end
            end

        end

        function [X_ECEF, t_k, A_k, e_n, E_k, u_k, i_k, OMEGA_k, PHI_k, r_k] ...
                = eph2xyz(obj, eph, ttx)
            import org.orekit.utils.*;
            mu = Constants.WGS84_EARTH_MU;
            omega_e = 0;%Constants.WGS84_EARTH_ANGULAR_VELOCITY; % Set to zero in this context because we are working in ECI
            % instead of ECEF coordinates. To work in ECEF, add omega_e to
            % the list of globals.

            % -------------------------------------------------------------------------
            % EXTRACT KEPLERIAN ORBITAL ELEMENTS
            % -------------------------------------------------------------------------

            % Find t_k.
            t_k = ttx - eph.Toe;

            % Account for week crossovers. Not needed in this analysis.
            % if t_k  > 302400
            %     t_k = t_k - 604800;
            % elseif t_k < -302400
            %     t_k = t_k + 604800;
            % end

            % Semi-major axis at reference time.
            A_0 = eph.Asqrt^2; % [m]
            A_k = A_0;

            % Mean motion.
            n_0 = sqrt(mu / A_0^3); % [rad/s]
            n_a = n_0 + eph.Delta_n;

            % Eccentricity.
            e_n = eph.e;

            % Corrected mean anomaly.
            M_k = eph.M0 + n_a*t_k;

            % Eccentric anomaly.
            E_k = obj.Keplers_Eqn(M_k,e_n);

            % True anomaly.
            tan_a = sqrt(1-e_n^2)*sin(E_k)/(1-e_n*cos(E_k));
            tan_b = (cos(E_k)-e_n)/(1-e_n*cos(E_k));

            v_k = atan2(tan_a,tan_b);

            % Harmonic correction terms.
            PHI_k = v_k + eph.Omega;

            du_k = eph.Cus*sin(2*PHI_k) + eph.Cuc*cos(2*PHI_k);
            dr_k = eph.Crs*sin(2*PHI_k) + eph.Crc*cos(2*PHI_k);
            di_k = eph.Cis*sin(2*PHI_k) + eph.Cic*cos(2*PHI_k);

            % Argument of latitude.
            u_k = PHI_k + du_k;

            % Radius.
            r_k = A_k*(1-e_n*cos(E_k)) + dr_k;

            % Inclination.
            i_k = eph.i0 + eph.IDOT*t_k + di_k;

            % Longitude of the ascending node.
            OMEGA_dot = eph.Omega_dot;
            OMEGA_k = eph.Omega0 + (OMEGA_dot - omega_e)*t_k  - omega_e*eph.Toe;


            % -------------------------------------------------------------------------
            % POSITION OF SATELLITE
            % -------------------------------------------------------------------------

            % SV position in the satellite orbit.
            xk_p = r_k*cos(u_k);
            yk_p = r_k*sin(u_k);

            % SV position in the ECEF frame.
            xk = xk_p*cos(OMEGA_k) - yk_p*cos(i_k)*sin(OMEGA_k);
            yk = xk_p*sin(OMEGA_k) + yk_p*cos(i_k)*cos(OMEGA_k);
            zk = yk_p*sin(i_k);

            X_ECEF = [xk yk zk];
            % -------------------------------------------------------------------------
            % CLOCK CORRECTION
            % -------------------------------------------------------------------------
            % Not used here.
            % t = ttx - eph.Toc;
            % dtsve = eph.a_f0 + eph.a_f1*t + eph.a_f2*t^2 + F*e_n*eph.Asqrt*sin(E_k);
        end

        function [a, ecc, inc, RAAN, omega, M0,...
                Cus, Cuc, Crs, Crc, Cis, Cic, ...
                IDOT, OMEGA_DOT, delta_n, flag,NumIter] = ...
                COE15_estimator(obj, time, pos, initial_guess, Wmat, ConvCrit, ...
                fit_parameters)
            % import org.orekit.utils.*;
            % R_e = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
            % mu = Constants.WGS84_EARTH_MU;
            % omega_e = Constants.WGS84_EARTH_ANGULAR_VELOCITY;

            m = length(time);
            % Rearrange the data - stack position vectors on top of each other.
            count = 1;
            for i = 1:m
                X_data_rearrange(count:count+2) = pos(i,:)';
                count = count + 3;
            end
            % Extract initial guess info.
            a     = initial_guess(1);    % [m]
            ecc   = initial_guess(2);    % [-]
            inc   = initial_guess(3);    % [rad]
            RAAN  = initial_guess(4);    % [rad]
            omega = initial_guess(5);    % [rad]
            M0    = initial_guess(6);    % [rad]
            Cus   = initial_guess(7);    % [rad]
            Cuc   = initial_guess(8);    % [rad]
            Crs   = initial_guess(9);    % [m]
            Crc   = initial_guess(10);   % [m]
            Cis   = initial_guess(11);   % [rad]
            Cic   = initial_guess(12);   % [rad]
            IDOT      = initial_guess(13); % [rad/s]
            OMEGA_DOT = initial_guess(14); % [rad/s]
            delta_n   = initial_guess(15); % [rad/s]
            % Scale to normalize lengths.
            scale_meters = a;
            % Non-dimensionalize a for numerical stability.
            a = a / scale_meters;
            Crs = Crs / scale_meters;
            Crc = Crc / scale_meters;
            % Form intial guess state vector.
            COEvec = [a, ecc, inc, RAAN, omega, M0, ...
                Cus, Cuc, Crs, Crc, Cis, Cic, ...
                IDOT, OMEGA_DOT, delta_n];
            % COEvec = initial_guess;
            % Logic variables.
            true  = 1;
            false = 0;
            done  = false;
            % Keep track of the number of iterations.
            NumIter = 0;
            % Set the failure flag to null. This will ultimately tell us whether or not
            % we have failed in producing a message.
            flag = 0;
            % Form weighting matrix based on diagonal input weights.
            count = 1;
            W = zeros(length(Wmat)*3);
            for i = 1:length(Wmat)
                W(count:count+2,count:count+2) = eye(3)*Wmat(i);
                count = count + 3;
            end
            % set the maximum number of iterations
            MaxIter = 100;
            % Define the error smoothing parameters.
            forward_diff_coeff = [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6];
            D_hat = zeros( 3 * m );
            num_neighbours = 6;
            D_hat = D_hat + diag( ones( 3 * m, 1 ) ) * forward_diff_coeff(1) ;
            for i = 1:num_neighbours
                D_hat = D_hat + diag( ones( 3 * (m - i), 1 ), 3 * i ) * ...
                    forward_diff_coeff(i+1);
            end
            D_hat = D_hat(1:end-3*num_neighbours, 1:end);
            u_hat = 0;
            % START ITERATION SCHEME FOR 15 ELEMENT ESTIMATE.
            while done == false && NumIter <= MaxIter
                % Initialize count variable.
                count = 1;
                for i = 1:length(time)
                    % Make an ephemeris structure to pass to the eph2xyz function.
                    eph.Asqrt  = sqrt(a * scale_meters);
                    eph.e      = ecc;
                    eph.i0     = inc;
                    eph.Omega0 = RAAN; % TODO this may need to be adjusted???
                    eph.Omega  = omega;
                    eph.M0     = M0;

                    eph.Cus    = Cus;
                    eph.Cuc    = Cuc;
                    eph.Crs    = Crs * scale_meters;
                    eph.Crc    = Crc * scale_meters;
                    eph.Cis    = Cis;
                    eph.Cic    = Cic;
                    eph.IDOT      = IDOT;
                    eph.Omega_dot = OMEGA_DOT;
                    eph.Delta_n   = delta_n;

                    eph.Toe = time(1);

                    % Define the time of transmission.
                    ttx = time(i);

                    % Compute the ECEF position as well as the other needed parameters.
                    [X_ECEF, t_k, A_k, e_n, E_k, u_k, i_k, OMEGA_k, PHI_k, r_k] = ...
                        obj.eph2xyz(eph, ttx);
                    A_k = A_k / scale_meters;
                    % Assign the theoretical position vector.
                    X_theo(count:count+2) = X_ECEF' / scale_meters;
                    % Compute vectors for partial derivatives.
                    r_k_hat = X_ECEF' / r_k;
                    r_k = r_k / scale_meters;
                    dr_k_hat_du_k = [
                        -sin(u_k) * cos(OMEGA_k) - cos(u_k) * cos(i_k) * sin(OMEGA_k);
                        -sin(u_k) * sin(OMEGA_k) + cos(u_k) * cos(i_k) * cos(OMEGA_k);
                        cos(u_k) * sin(i_k)];
                    dr_k_hat_du_k2 = [
                        -sin(u_k) * sin(OMEGA_k) - cos(u_k) * cos(i_k) * sin(OMEGA_k);
                        -sin(u_k) * sin(OMEGA_k) + cos(u_k) * cos(i_k) * cos(OMEGA_k);
                        cos(u_k) * sin(i_k)];
                    dr_k_hat_du_k2 = dr_k_hat_du_k;
                    dr_k_hat_di_k = sin(u_k) * [
                        sin(i_k) * sin(OMEGA_k);
                        -sin(i_k) * cos(OMEGA_k);
                        cos(i_k)];

                    % Compute derivatives.
                    drda = (1 - e_n * cos(E_k)) * r_k_hat;
                    drde = A_k * (e_n - cos(E_k)) / (1 - e_n * cos(E_k)) * r_k_hat +...
                        r_k * (2 * sin(E_k) - e_n * sin(E_k) * cos(E_k) - e_n ^ 2 * sin(E_k)) / ...
                        sqrt(1 - e_n ^ 2) / (1 - e_n * cos(E_k)) ^ 2 * dr_k_hat_du_k;
                    drdi0 = r_k * dr_k_hat_di_k;
                    drdOMEGA0 = r_k * [
                        -cos(u_k) * sin(OMEGA_k) - sin(u_k) * cos(i_k) * cos(OMEGA_k);
                        cos(u_k) * cos(OMEGA_k) - sin(u_k) * cos(i_k) * sin(OMEGA_k);
                        0];
                    drdomega = r_k * dr_k_hat_du_k2;
                    drdM0 = A_k * e_n * sin(E_k) / (1 - e_n * cos(E_k)) * r_k_hat + ...
                        r_k * sqrt(1 - e_n^2) / (1 - e_n * cos(E_k)) ^2 * dr_k_hat_du_k;
                    drdCus = sin(2*PHI_k) * r_k * dr_k_hat_du_k2;
                    drdCuc = cos(2*PHI_k) * r_k * dr_k_hat_du_k2;
                    drdCrs = sin(2*PHI_k) * r_k_hat;
                    drdCrc = cos(2*PHI_k) * r_k_hat;
                    drdCis = sin(2*PHI_k) * r_k * dr_k_hat_di_k;
                    drdCic = cos(2*PHI_k) * r_k * dr_k_hat_di_k;
                    drdIDOT = t_k * drdi0;
                    drdOMEGA_DOT = t_k * drdOMEGA0;
                    drddelta_n = t_k * drdM0;
                    % Form Jacobian matrix.
                    A(count:count+2,:) = [
                        drda, drde, drdi0, drdOMEGA0, drdomega, drdM0, ...
                        drdCus, drdCuc, drdCrs, drdCrc, drdCis, drdCic, ...
                        drdIDOT, drdOMEGA_DOT, drddelta_n].*...
                        [fit_parameters;fit_parameters;fit_parameters];
                    % Update the counter.
                    count = count + 3;
                end

                % Solve for the update using Matlab's matrix divide.
                y = X_data_rearrange' / scale_meters - X_theo';
                dCOEvec = A \ y;
                % If we want weighted least squares.
                dCOEvec = (W * A) \ ( W * y);
                % Update orbital element vector.
                COEvec = COEvec + dCOEvec';
                % Assign the new values of the orbital elements.
                a     = COEvec(1);    % [m]
                ecc   = COEvec(2);    % [-]
                inc   = COEvec(3);    % [rad]
                RAAN  = COEvec(4);    % [rad]
                omega = COEvec(5);    % [rad]
                M0    = COEvec(6);    % [rad]
                Cus   = COEvec(7);    % [rad]
                Cuc   = COEvec(8);    % [rad]
                Crs   = COEvec(9);    % [m]
                Crc   = COEvec(10);   % [m]
                Cis   = COEvec(11);   % [rad]
                Cic   = COEvec(12);   % [rad]
                IDOT      = COEvec(13); % [rad/s]
                OMEGA_DOT = COEvec(14); % [rad/s]
                delta_n   = COEvec(15); % [rad/s]
                % Mitigate negative eccetricity.
                if ecc < 0
                    ecc = abs( ecc );
                    COEvec(2) = ecc;
                    M0 = M0 + pi;
                    COEvec(6) = M0;
                    omega = omega + pi;
                    COEvec(5) = omega;
                end
                % If we get nonsensical things, we have failed.
                if a < 0 || a > 1.5 || ecc > 1
                    done = true;
                    NumIter = MaxIter;
                    flag = 1;
                end
                % Update number of iterations.
                NumIter = NumIter + 1;

                % Check for convergence.
                newton_decrement = norm( A * dCOEvec );
                if newton_decrement < ConvCrit
                    done = true;
                end
            end
            % Re-dimensionalize.
            a = a * scale_meters; % [m]
            Crs = Crs * scale_meters; % [m]
            Crc = Crc * scale_meters; % [m]
            % Put the angular quantities in the correct range (between -pi and pi)
            inc = wrapToPi(inc);
            RAAN = wrapToPi(RAAN);
            omega = wrapToPi(omega);
            M0 = wrapToPi(M0);
            % Final check for failure.
            if NumIter >= MaxIter
                % Set flag.
                flag = 1;
                % Set output to NaN.
                a     = NaN;
                ecc   = NaN;
                inc   = NaN;
                RAAN  = NaN;
                omega = NaN;
                M0    = NaN;
                Cus   = NaN;
                Cuc   = NaN;
                Crs   = NaN;
                Crc   = NaN;
                Cis   = NaN;
                Cic   = NaN;
                IDOT  = NaN;
                OMEGA_DOT = NaN;
                delta_n   = NaN;
                NumIter   = NaN;
            end
        end

        function [a, ecc, inc, RAAN, omega, M0,...
                Cus, Cuc, Crs, Crc, Cis, Cic, ...
                IDOT, OMEGA_DOT, delta_n, flag, NumIter, fit_type] = ...
                COE15_estimator_wrapper(obj, time, pos, vel,...
                initial_guess, Wmat, ConvCrit, ...
                fit_parameters, theta_g, coeff_R2, coeff_A2, coeff_C2)
            % Define the fit parameter scenarios to be tested.
            fit_scenarios = repmat(fit_parameters, 10, 1);
            % No RAAN.
            fit_scenarios(2, 4) = 0;
            % No omega.
            fit_scenarios(3, 5) = 0;
            % No inclination.
            fit_scenarios(4, 3) = 0;
            % No RAAN / inc.
            fit_scenarios(5, 3) = 0;
            fit_scenarios(5, 4) = 0;
            % No omega / inc.
            fit_scenarios(6, 3) = 0;
            fit_scenarios(6, 5) = 0;
            % No RAAN / ecc.
            fit_scenarios(7, 2) = 0;
            fit_scenarios(7, 4) = 0;
            % No omega / ecc.
            fit_scenarios(8, 2) = 0;
            fit_scenarios(8, 4) = 0;
            % No RAAN / omega.
            fit_scenarios(9, 4) = 0;
            fit_scenarios(9, 5) = 0;
            % No RAAN / omega / inc.
            fit_scenarios(10, 3) = 0;
            fit_scenarios(10, 4) = 0;
            fit_scenarios(10, 5) = 0;
            % No RAAN / omega / ecc.
            fit_scenarios(11, 2) = 0;
            fit_scenarios(11, 4) = 0;
            fit_scenarios(11, 5) = 0;
            % No RAAN / omega / ecc / inc.
            fit_scenarios(12, 2) = 0;
            fit_scenarios(12, 3) = 0;
            fit_scenarios(12, 4) = 0;
            fit_scenarios(12, 5) = 0;
            % Get the number of scenarios.
            [num_scenarios,~] = size(fit_scenarios);
            % Initialize the rms ure (this is a dummy value).
            rms_ure_best = 999;
            for i = 1:num_scenarios
                % Try fitting without removing elements.
                [a_test, ecc_test, inc_test, RAAN_test, omega_test, M0_test,...
                    Cus_test, Cuc_test, Crs_test, Crc_test, Cis_test, Cic_test, ...
                    IDOT_test, OMEGA_DOT_test, delta_n_test, ...
                    flag_test, NumIter_test] = ...
                    obj.COE15_estimator(time, pos, initial_guess, Wmat, ConvCrit, ...
                    fit_scenarios(i,:));

                % If the inclination is less than zero, we have failed.
                if inc_test < 0
                    flag_test = 1;
                end

                % Evaluate the message performance.
                if flag_test == 0
                    [~, ~, rms_ure, ~, ~, ~, ~] = ...
                        obj.eph_error_analysis(sqrt(a_test), ecc_test, inc_test,...
                        RAAN_test, omega_test, M0_test, ...
                        Cus_test, Cuc_test, Crc_test, Crs_test, ...
                        Cic_test, Cis_test, ...
                        IDOT_test, OMEGA_DOT_test, delta_n_test, ...
                        time, pos, vel, theta_g, ...
                        coeff_R2, coeff_A2, coeff_C2);
                end

                if flag_test == 0
                    if rms_ure < rms_ure_best
                        % Assign the best rms_ure seen so far.
                        rms_ure_best = rms_ure;

                        % Assign the fit type.
                        fit_type = i;

                        % Assign the orbital elements as the output.
                        a = a_test; % [m]
                        ecc = ecc_test; % [-]
                        inc = inc_test; % [rad]
                        RAAN = RAAN_test; % [rad]
                        omega = omega_test; % [rad]
                        M0 = M0_test; % [rad]

                        Cus   = Cus_test;    % [rad]
                        Cuc   = Cuc_test;    % [rad]
                        Crs   = Crs_test;    % [m]
                        Crc   = Crc_test;   % [m]
                        Cis   = Cis_test;   % [rad]
                        Cic   = Cic_test;   % [rad]

                        IDOT      = IDOT_test; % [rad/s]
                        OMEGA_DOT = OMEGA_DOT_test; % [rad/s]
                        delta_n   = delta_n_test; % [rad/s]

                        % Number of iterations.
                        NumIter = NumIter_test;

                        % Set the flag.
                        flag = flag_test;

                        % If we succeed without having to remove any elements, stop.
                        if i == 0
                            break;
                        end
                    end
                end
            end
        end

        function [error_3d, rms_error, rms_ure, ...
                error_radial, error_along_track, error_cross_track, ...
                eph] = eph_error_analysis(obj, sqrt_a, ecc, inc, RAAN, omega, M0, ...
                Cus, Cuc, Crc, Crs, Cic, Cis, ...
                IDOT, OMEGA_DOT, delta_n, time, pos, vel, theta_g, ...
                coeff_R2, coeff_A2, coeff_C2)
            import org.orekit.utils.*;
            % global omega_e
            omega_e = 0;%Constants.WGS84_EARTH_ANGULAR_VELOCITY; % Set to zero in this context because we are working in ECI
            % instead of ECEF coordinates. To work in ECEF, add omega_e to
            % the list of globals.
            % Make an ephemeris structure to pass to the eph2xyz function.
            eph.Asqrt  = sqrt_a;
            eph.e      = ecc;
            eph.i0     = inc;
            eph.Omega0 = RAAN;
            eph.Omega  = omega;
            eph.M0     = M0;
            eph.Cus    = Cus;
            eph.Cuc    = Cuc;
            eph.Crs    = Crs;
            eph.Crc    = Crc;
            eph.Cis    = Cis;
            eph.Cic    = Cic;
            eph.IDOT      = IDOT;
            eph.Omega_dot = OMEGA_DOT;
            eph.Delta_n   = delta_n;
            eph.Toe = time(1);
            % Number of data points.
            n_data = length(time);
            % Initialize vectors.
            pos_theo = NaN(n_data, 3);
            error_3d = NaN(n_data, 1);
            error_radial = NaN(n_data, 1);
            error_along_track = NaN(n_data, 1);
            error_cross_track = NaN(n_data, 1);
            % Compute message error.
            for i = 1:n_data
                % Transmission time.
                ttx = time(i);
                % Compute message theorectical position.
                [x_ecef_eph, ~, ~, ~, ~, ~, ~, ~, ~, ~] ...
                    = obj.eph2xyz(eph, ttx);
                pos_theo(i,:) = x_ecef_eph;
                % Get the sidereal time.
                theta_g_epoch = theta_g + omega_e * time(i);
                % Convert ECEF pos / vel to ECI
                [x_eci, v_eci] = obj.ECEF2ECI( pos(i,:)', vel(i,:)', theta_g_epoch );
                [x_eci_eph, ~] = obj.ECEF2ECI( x_ecef_eph', NaN(3,1), theta_g_epoch );
                % Compute transformation matrix ECI_2_RIC (radial / along track /
                % cross track).
                ECI_2_RIC_Mat = obj.ECI2RIC( x_eci, v_eci );
                x_ric = ECI_2_RIC_Mat * x_eci;
                x_ric_eph = ECI_2_RIC_Mat * x_eci_eph;
                % Compute the radial / along-track / cross-track errors.
                error_radial(i) = x_ric_eph(1) - x_ric(1);
                error_along_track(i) = x_ric_eph(2) - x_ric(2);
                error_cross_track(i) = x_ric_eph(3) - x_ric(3);
                % Compute 3D error.
                error_3d(i) = norm(x_ecef_eph - pos(i,:));

            end

            % RMS 3D error.
            rms_error = ( sum( ...
                error_radial.^2 + error_along_track.^2 + error_cross_track.^2 ...
                ) / n_data ) ^ (1/2);

            % Calculate the RMS error component errors.
            RMS_radial = rms( error_radial );
            RMS_cross_track = rms( error_cross_track );
            RMS_along_track = rms( error_along_track );

            % Calculate the RMS orbit only ure.
            rms_ure = sqrt( coeff_A2 * RMS_along_track ^ 2 + ...
                coeff_C2 * RMS_cross_track ^ 2 + ...
                coeff_R2 * RMS_radial ^ 2 );
        end

        function lgnssEphemeris = getNavigationMessageList(obj, ...
                thisOrbit, wholeDurationInhours, eachMessageFitIntervals, timeBetweenMessages)
            import org.orekit.utils.*
            if wholeDurationInhours * 3600.0 < eachMessageFitIntervals
                return
            end
            R_e = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
            mu = Constants.WGS84_EARTH_MU;
            midProp = obj.getPropagator("Get", thisOrbit, "H");
            ephGenerator = midProp.getEphemerisGenerator();
            midProp.propagate(obj.startTime, obj.startTime.shiftedBy(wholeDurationInhours*3600.0));
            thisPropagator = ephGenerator.getGeneratedEphemeris();
            currentDateTime = obj.startTime;
            trueEphemeris = zeros(floor((wholeDurationInhours*3600)/obj.plotStepTime) + 1,7);
            mn = 0;
            while currentDateTime.compareTo(obj.startTime.shiftedBy(wholeDurationInhours*3600.0)) <=0
                mn = mn + 1;
                trueEphemeris(mn , 1) = mn - 1;
                truePV = thisPropagator.propagate(currentDateTime).getPVCoordinates();
                trueEphemeris(mn , 2 : 7) = [truePV.getPosition().getX(), truePV.getPosition().getY(), ...
                    truePV.getPosition().getZ(), truePV.getVelocity().getX(), ...
                    truePV.getVelocity().getY(), truePV.getVelocity().getZ()];
                currentDateTime = currentDateTime.shiftedBy(obj.plotStepTime);
            end
            if (mn) < size(trueEphemeris , 1)
                ...
            end
            orbit_data.pos_m = trueEphemeris(:,2:4); %ECI
            orbit_data.vel_m_s = trueEphemeris(:,5:7); %ECI
            orbit_data.elapsed_time_sec = trueEphemeris(:,1);
            ICUIS = obj.ECI2COE(trueEphemeris(1,2:4),trueEphemeris(1,5:7));
            SMA = ICUIS.a; % [m]
            ecc = ICUIS.e; % [-]
            inc = ICUIS.i; % [deg]
            RAAN = ICUIS.RAAN; % [deg]
            AOP = ICUIS.omega; % [deg]
            MA = ICUIS.M; % [deg]

            % Fit interval of interest.
            fit_interval = eachMessageFitIntervals; % [seconds]

            % File names.
            file_altitudes = SMA - R_e; % [km]

            % Total time in file.
            time_in_file = orbit_data.elapsed_time_sec(end); % [sec]

            % Start the clock.
            tic;

            % Get the analytical coefficients for the URE equations.
            [coeff_A2, coeff_C2, coeff_R2, coeff_T2, coeff_RT, theta] = ...
                obj.analytic_URE_eqn(  file_altitudes );

            % Determine the message start times.
            message_start_times = ...
                0:timeBetweenMessages:(time_in_file - fit_interval);

            % Determine the number of messages that can be made per file.
            num_eph_per_file = length( message_start_times );

            % Initialize variables.
            rms_ure_Save = NaN(num_eph_per_file, 1);
            rms_3D_Save = NaN(num_eph_per_file, 1);
            max_3D_Save = NaN(num_eph_per_file, 1);
            Num_Iter_Save = NaN(num_eph_per_file, 1);
            convergence_crit_Save = NaN(num_eph_per_file, 1);
            failure_flag_Save = zeros(num_eph_per_file, 1);
            fit_type_Save = zeros(num_eph_per_file, 1);
            eph_datenum_Save = NaN(num_eph_per_file, 1);

            eph_Save(num_eph_per_file).Asqrt  = [];
            eph_Save(num_eph_per_file).e      =[];
            eph_Save(num_eph_per_file).i0     =[];
            eph_Save(num_eph_per_file).Omega0 =[];
            eph_Save(num_eph_per_file).Omega  =[];
            eph_Save(num_eph_per_file).M0     =[];

            eph_Save(num_eph_per_file).Cus    =[];
            eph_Save(num_eph_per_file).Cuc    =[];
            eph_Save(num_eph_per_file).Crs    =[];
            eph_Save(num_eph_per_file).Crc    =[];
            eph_Save(num_eph_per_file).Cis    =[];
            eph_Save(num_eph_per_file).Cic    =[];

            eph_Save(num_eph_per_file).IDOT      =[];
            eph_Save(num_eph_per_file).Omega_dot =[];
            eph_Save(num_eph_per_file).Delta_n   =[];
            eph_Save(num_eph_per_file).Toe   =[];

            % Start threshold for progress bar.
            progress = 0; % [percent]

            for idx_message = 1:num_eph_per_file
                % Percent done.
                percent_done = floor(idx_message/num_eph_per_file * 100);
                if percent_done > progress
                    disp(['Percent done: ', num2str(percent_done), ', Elapsed time is: ', num2str(toc/60), ' [min]']);
                    progress = percent_done;
                end

                % Define the time vector for fitting.
                time = ...
                    orbit_data.elapsed_time_sec(...
                    orbit_data.elapsed_time_sec <= fit_interval );

                % Get the start / end index for the data to fit to.
                idx_start = find(orbit_data.elapsed_time_sec == ...
                    message_start_times(idx_message));
                idx_end = find(orbit_data.elapsed_time_sec == ...
                    message_start_times(idx_message) + fit_interval);

                % Get position and velocity vectors
                pos = ...
                    orbit_data.pos_m(idx_start:idx_end,:);
                vel = ...
                    orbit_data.vel_m_s(idx_start:idx_end,:);

                % Convert ECEF to ECI coordinates.
                % NOTE: We'll work in ECI coordinates for the purposes of this
                %       experiement but, ECEF is also possible here with some
                %       small changes.
                % theta_g = utc2gmst( datevec(orbit_data.datenum(idx_start)) ); % [rad]
                % [R_test, V_test] = ECEF2ECI(...
                %     pos(1,:)',  vel(1,:)', theta_g)

                % Since we're dealing with ECI vs ECEF coordinates, we'll use a
                % zero offset between them (for the purposes of reusing other
                % code).
                theta_g = 0; % [rad]

                % Form initial guess with the 6 Keplerian elements.
                [coe, ~, ~] = obj.ECI2COE(pos(1,:), vel(1,:));

                % Form the initial guess for the estimator.
                a = coe.a; % [m]
                n = sqrt(mu/a^3); % [rad/sec]
                ecc = coe.e; % [-]
                inc = coe.i * pi / 180; % [rad]
                RAAN = coe.RAAN * pi / 180; % [rad]
                omega = coe.omega * pi /180; % [rad]
                M0 = coe.M * pi / 180; % [rad]

                Cus = 0; % [rad]
                Cuc = 0; % [rad]
                Crs = 0; % [rad]
                Crc = 0; % [rad]
                Cis = 0; % [rad]
                Cic = 0; % [rad]

                IDOT = 0; % [rad/s]
                OMEGA_DOT = 0; % [rad/s]
                delta_n = 0; % [rad/s]

                % Form the initial guess.
                initial_guess(1) = a;
                initial_guess(2) = ecc;
                initial_guess(3) = inc;
                initial_guess(4) = RAAN;
                initial_guess(5) = omega;
                initial_guess(6) = M0;

                initial_guess(7)  = Cus;
                initial_guess(8)  = Cuc;
                initial_guess(9)  = Crs;
                initial_guess(10) = Crc;
                initial_guess(11) = Cis;
                initial_guess(12) = Cic;

                initial_guess(13) = IDOT;
                initial_guess(14) = OMEGA_DOT;
                initial_guess(15) = delta_n;

                % Define the convergence criteria.
                ConvCrit = 1e-11;

                % Fit ephemeris parameters or subset.
                % Define the weighting matrix, use the identity matrix for now.
                Wmat = ones( size(time) );

                fit_parameters = zeros(1,15);

                % Keplerian Elements.
                fit_parameters(1) = 1; % a
                fit_parameters(2) = 1; % e
                fit_parameters(3) = 1; % inc
                fit_parameters(4) = 1; % RAAN
                fit_parameters(5) = 1; % omega
                fit_parameters(6) = 1; % M0

                % Corrections.
                fit_parameters(7) = 1; % Cus
                fit_parameters(8) = 1; % Cuc
                fit_parameters(9) = 1; % Crs
                fit_parameters(10) = 1; % Crc
                fit_parameters(11) = 1; % Cis
                fit_parameters(12) = 1; % Cic
                fit_parameters(13) = 1; % IDOT
                fit_parameters(14) = 1; % OMEGA_DOT
                fit_parameters(15) = 1; % delta_n

                % Fit parameters.
                [a, ecc, inc, RAAN, omega, M0,...
                    Cus, Cuc, Crs, Crc, Cis, Cic, ...
                    IDOT, OMEGA_DOT, delta_n, flag, NumIter, fit_type] = ...
                    obj.COE15_estimator_wrapper(time, pos, vel, initial_guess, ...
                    Wmat, ConvCrit, fit_parameters, ...
                    theta_g, coeff_R2, coeff_A2, coeff_C2);

                % Quantize message parameters.
                %     [a, ecc, inc, RAAN, omega, M0, Cus, Cuc, IDOT, numbits] = ...
                %         bit_reduction(a, ecc, inc, RAAN, omega, M0, Cus, Cuc, IDOT);

                % Error analysis.
                [error_3d, rms_error, rms_ure, ...
                    error_radial, error_along_track, error_cross_track, ...
                    eph] = obj.eph_error_analysis(sqrt(a), ecc, inc, RAAN, omega, M0, ...
                    Cus, Cuc, Crc, Crs, Cic, Cis, ...
                    IDOT, OMEGA_DOT, delta_n, time, pos, vel, theta_g, ...
                    coeff_R2, coeff_A2, coeff_C2);
                % Save results.
                rms_ure_Save(idx_message) = rms_ure;
                rms_3D_Save(idx_message) = rms_error;
                max_3D_Save(idx_message) = max(error_3d);
                Num_Iter_Save(idx_message) = NumIter;
                convergence_crit_Save(idx_message) =  ...
                    ConvCrit;
                failure_flag_Save(idx_message) = flag;
                fit_type_Save(idx_message) = fit_type;
                eph_Save(idx_message) = eph;
            end % end idx_message
            lgnssEphemeris.rms_ure_Save = rms_ure_Save;
            lgnssEphemeris.rms_3D_Save = rms_3D_Save;
            lgnssEphemeris.max_3D_Save = max_3D_Save;
            lgnssEphemeris.Num_Iter_Save = Num_Iter_Save;
            lgnssEphemeris.convergence_crit_Save = convergence_crit_Save;
            lgnssEphemeris.failure_flag_Save = failure_flag_Save;
            lgnssEphemeris.fit_type_Save = fit_type_Save;
            lgnssEphemeris.eph_Save = eph_Save;
        end

        function ephemerisError = getLgnssError(obj, ephemerisStruct, trueOrbit, wholeDurationInhours, errorStep)

            initialDate = trueOrbit.getDate();
            if initialDate.durationFrom(ephemerisStruct.Epoch) ~= 0
                ephemerisError = [];
                return
            end


            midProp = obj.getPropagator("Get", trueOrbit, "H");
            ephGenerator = midProp.getEphemerisGenerator();
            midProp.propagate(initialDate, initialDate.shiftedBy(wholeDurationInhours*3600.0));
            thisPropagator = ephGenerator.getGeneratedEphemeris();
            currentDateTime = initialDate;
            ephemerisError = zeros(floor((wholeDurationInhours*3600)/errorStep) + 1,4);
            mn = 0;
            while currentDateTime.compareTo(initialDate.shiftedBy(wholeDurationInhours*3600.0)) <=0
                mn = mn + 1;
                truePV = thisPropagator.propagate(currentDateTime).getPVCoordinates();
                trueEphemeris = [truePV.getPosition().getX(), truePV.getPosition().getY(), ...
                    truePV.getPosition().getZ()];
                ephemerisError(mn , 1 : 3) = trueEphemeris - obj.eph2xyz(ephemerisStruct, (mn-1)*errorStep );
                ephemerisError(mn , 4) = norm(ephemerisError(mn , 1 : 3));
                currentDateTime = currentDateTime.shiftedBy(errorStep);
            end

        end

        function lgnssEphemeris = getNavigationMessage(obj, ...
                thisOrbit, messageFitInterval)
            initialDate = thisOrbit.getDate();
            observationTimeStep = 1.0;
            import org.orekit.utils.*
            R_e = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
            mu = Constants.WGS84_EARTH_MU;
            midProp = obj.getPropagator("Get", thisOrbit, "H");
            ephGenerator = midProp.getEphemerisGenerator();
            midProp.propagate(initialDate, initialDate.shiftedBy(messageFitInterval));
            thisPropagator = ephGenerator.getGeneratedEphemeris();
            trueEphemeris = zeros(floor((messageFitInterval/observationTimeStep)) + 1,7);
            trueEphemerisEcef = zeros(floor((messageFitInterval/observationTimeStep)) + 1,7);
            unixTimes = zeros(floor((messageFitInterval/observationTimeStep)) + 1,1);
            currentDateTime = initialDate;
            mn = 0;
            while currentDateTime.compareTo(initialDate.shiftedBy(messageFitInterval)) <=0
                mn = mn + 1;
                trueEphemeris(mn , 1) = mn - 1;
                trueEphemerisEcef(mn , 1) = mn - 1;
                unixTimes(mn , 1) = obj.getUnixTime(currentDateTime);
                truePV = thisPropagator.propagate(currentDateTime).getPVCoordinates();
                trueEphemeris(mn , 2 : 7) = [truePV.getPosition().getX(), truePV.getPosition().getY(), ...
                    truePV.getPosition().getZ(), truePV.getVelocity().getX(), ...
                    truePV.getVelocity().getY(), truePV.getVelocity().getZ()];

                truePVEcef = thisPropagator.propagate(currentDateTime).getPVCoordinates(obj.ecefFrame);
                trueEphemerisEcef(mn , 2 : 7) = [truePVEcef.getPosition().getX(), truePVEcef.getPosition().getY(), ...
                    truePVEcef.getPosition().getZ(), truePVEcef.getVelocity().getX(), ...
                    truePVEcef.getVelocity().getY(), truePVEcef.getVelocity().getZ()];
                currentDateTime = currentDateTime.shiftedBy(observationTimeStep);
            end
            if (mn) < size(trueEphemeris , 1)
                ...
            end
            orbit_data.pos_m = trueEphemeris(:,2:4); %ECI
            orbit_data.vel_m_s = trueEphemeris(:,5:7); %ECI
            orbit_data.elapsed_time_sec = trueEphemeris(:,1);
            ICUIS = obj.ECI2COE(trueEphemeris(1,2:4),trueEphemeris(1,5:7));
            SMA = ICUIS.a; % [m]

            % Fit interval of interest.
            fit_interval = messageFitInterval; % [seconds]

            % File names.
            file_altitudes = SMA - R_e; % [km]

            % Total time in file.
            time_in_file = orbit_data.elapsed_time_sec(end); % [sec]

            % Start the clock.
            tic;

            % Get the analytical coefficients for the URE equations.
            [coeff_A2, coeff_C2, coeff_R2, ~, ~, ~] = ...
                obj.analytic_URE_eqn(  file_altitudes );

            % Determine the message start times.
            message_start_times = ...
                0:messageFitInterval:(time_in_file - fit_interval);

            % Determine the number of messages that can be made per file.
            num_eph_per_file = length( message_start_times );

            % Initialize variables.

            lgnssEphemeris(num_eph_per_file).Asqrt  = [];
            lgnssEphemeris(num_eph_per_file).e      =[];
            lgnssEphemeris(num_eph_per_file).i0     =[];
            lgnssEphemeris(num_eph_per_file).Omega0 =[];
            lgnssEphemeris(num_eph_per_file).Omega  =[];
            lgnssEphemeris(num_eph_per_file).M0     =[];

            lgnssEphemeris(num_eph_per_file).Cus    =[];
            lgnssEphemeris(num_eph_per_file).Cuc    =[];
            lgnssEphemeris(num_eph_per_file).Crs    =[];
            lgnssEphemeris(num_eph_per_file).Crc    =[];
            lgnssEphemeris(num_eph_per_file).Cis    =[];
            lgnssEphemeris(num_eph_per_file).Cic    =[];

            lgnssEphemeris(num_eph_per_file).IDOT      =[];
            lgnssEphemeris(num_eph_per_file).Omega_dot =[];
            lgnssEphemeris(num_eph_per_file).Delta_n   =[];
            lgnssEphemeris(num_eph_per_file).Toe   =[];
            lgnssEphemeris(num_eph_per_file).Epoch = [];

            idx_message = 1;

            % Define the time vector for fitting.
            time = ...
                orbit_data.elapsed_time_sec(...
                orbit_data.elapsed_time_sec <= fit_interval );

            % Get the start / end index for the data to fit to.
            idx_start = find(orbit_data.elapsed_time_sec == ...
                message_start_times(idx_message));
            idx_end = find(orbit_data.elapsed_time_sec == ...
                message_start_times(idx_message) + fit_interval);

            % Get position and velocity vectors
            pos = ...
                orbit_data.pos_m(idx_start:idx_end,:);
            vel = ...
                orbit_data.vel_m_s(idx_start:idx_end,:);

            theta_g = obj.utc2gmst( datevec(unixTimes(1)) ); % [rad]

            % Form initial guess with the 6 Keplerian elements.
            [coe, ~, ~] = obj.ECI2COE(trueEphemeris(1,2:4),trueEphemeris(1,5:7));

            % Form the initial guess for the estimator.
            a = coe.a; % [m]
            n = sqrt(mu/a^3); % [rad/sec]
            ecc = coe.e; % [-]
            inc = coe.i * pi / 180; % [rad]
            RAAN = coe.RAAN * pi / 180; % [rad]
            omega = coe.omega * pi /180; % [rad]
            M0 = coe.M * pi / 180; % [rad]

            Cus = 0; % [rad]
            Cuc = 0; % [rad]
            Crs = 0; % [rad]
            Crc = 0; % [rad]
            Cis = 0; % [rad]
            Cic = 0; % [rad]

            IDOT = 0; % [rad/s]
            OMEGA_DOT = 0; % [rad/s]
            delta_n = 0; % [rad/s]

            % Form the initial guess.
            initial_guess(1) = a;
            initial_guess(2) = ecc;
            initial_guess(3) = inc;
            initial_guess(4) = RAAN;
            initial_guess(5) = omega;
            initial_guess(6) = M0;

            initial_guess(7)  = Cus;
            initial_guess(8)  = Cuc;
            initial_guess(9)  = Crs;
            initial_guess(10) = Crc;
            initial_guess(11) = Cis;
            initial_guess(12) = Cic;

            initial_guess(13) = IDOT;
            initial_guess(14) = OMEGA_DOT;
            initial_guess(15) = delta_n;

            % Define the convergence criteria.
            ConvCrit = 1e-11;

            % Fit ephemeris parameters or subset.
            % Define the weighting matrix, use the identity matrix for now.
            Wmat = ones( size(time) );

            fit_parameters = ones(1,15);

            % Fit parameters.
            [a, ecc, inc, RAAN, omega, M0,...
                Cus, Cuc, Crs, Crc, Cis, Cic, ...
                IDOT, OMEGA_DOT, delta_n, ~, ~, ~] = ...
                obj.COE15_estimator_wrapper(time, pos, vel, initial_guess, ...
                Wmat, ConvCrit, fit_parameters, ...
                theta_g, coeff_R2, coeff_A2, coeff_C2);

            % Error analysis.
            [~, ~, ~, ...
                ~, ~, ~, ...
                lgnssEphemeris] = obj.eph_error_analysis(sqrt(a), ecc, inc, RAAN, omega, M0, ...
                Cus, Cuc, Crc, Crs, Cic, Cis, ...
                IDOT, OMEGA_DOT, delta_n, time, pos, vel, theta_g, ...
                coeff_R2, coeff_A2, coeff_C2);
            lgnssEphemeris(num_eph_per_file).Epoch = initialDate;
        end

        function linkedSatelliteStates = getLinkedSatellites(obj, pointLla, constellationOnEpoch, ...
                NoOfIntervals, wholeDurationInHours, maskInDegrees)
            intervalsDuration = (wholeDurationInHours*3600.0/NoOfIntervals);
            linkedSatelliteStates = cell(1,NoOfIntervals + 1);
            pointEcef = obj.llaToEcef(pointLla);
            currentDateTimes = constellationOnEpoch(1).getDate();
            propagators = obj.getConstellationPropagators(constellationOnEpoch, wholeDurationInHours*3600.0);
            for i = 1:NoOfIntervals + 1
                constellationOnEpoch = obj.getConstellaionStatesAt(propagators, ...
                    intervalsDuration*(i-1));
                constellationEcef = obj.getConstellationEcef(constellationOnEpoch);
                linkedSatelliteStates{1,i} = obj.getLinkedSatellitesEcef(constellationEcef,...
                    pointLla, maskInDegrees);
            end

        end

        function PsuedoRangeAndEcef = getPsuedoRangeAndEcefFromEph(obj, linkedSatelliteStates, EphemerisList,...
                pointEcef, NoOfIntervals, wholeDurationInHours)
            PsuedoRangeAndEcef = cell(1, NoOfIntervals + 1);
            intervalsDuration = (wholeDurationInHours*3600.0/(NoOfIntervals + 1));
            for i = 1:NoOfIntervals + 1
                linkedStates = linkedSatelliteStates{1,i};
                if isempty(linkedStates)
                else
                    indexArray = linkedStates(:,1);
                    states = zeros(size(indexArray,1), 4);
                    for j = 1:size(indexArray,1)
                        states(j,1:3) = obj.getEcefFromNavigationMessage(EphemerisList(indexArray(j)), intervalsDuration*(i-1));
                        states(j,4) = norm(states(j,1:3) - pointEcef);
                    end
                    PsuedoRangeAndEcef{1,i} = states;
                end
            end
        end

        function positionErrors = getUserErrorSingleMessage(obj, pointLla, constellationOnEpoch,...
                NoOfIntervals, wholeDurationInHours, maskInDegrees)
            intervalsDuration = (wholeDurationInHours*3600.0/NoOfIntervals);
            pointEcef = obj.llaToEcef(pointLla);
            if wholeDurationInHours*3600.0 > 30.0 * 60.0
                positionErrors = [];
                return
            end
            positionErrors = zeros(NoOfIntervals + 1, 4);
            lgnssEphemeris = [];
            for i = 1:size(constellationOnEpoch,1)
                lgnssEphemeris = [lgnssEphemeris; ...
                    obj.getNavigationMessage(constellationOnEpoch(i), wholeDurationInHours*3600.0)];
            end
            linkedSatelliteStates = obj.getLinkedSatellites(pointLla, ...
                constellationOnEpoch, NoOfIntervals, wholeDurationInHours, maskInDegrees);
            PsuedoRangeAndEcef = obj.getPsuedoRangeAndEcefFromEph(linkedSatelliteStates, lgnssEphemeris,...
                pointEcef, NoOfIntervals, wholeDurationInHours);
            for i = 1:NoOfIntervals + 1
                % linkedStates = linkedSatelliteStates{1,i};indexArray = linkedStates(:,1);
                linkedStates = PsuedoRangeAndEcef{1,i};
                if isempty(linkedStates)
                elseif size(linkedStates, 1) < 4
                else
                    positions = linkedStates(:,1:3)'; % u were here
                    ranges = linkedStates(:, 4)';
                    weights = diag(ones(1,length(ranges)));
                    [N1,~] = obj.trilateration(positions, ranges, weights);
                    estimatedPosition = N1(2:4,1)';
                    positionErrors(i,1:3) = estimatedPosition - pointEcef;
                    positionErrors(i,4) = norm(estimatedPosition - pointEcef);
                end
            end

        end

        function positionECEF = getEcefFromNavigationMessage(obj, ephemeris, interval)
            import org.hipparchus.geometry.euclidean.threed.*;
            posEci = obj.eph2xyz(ephemeris, interval);
            thisPositionECI = Vector3D(posEci(1),posEci(2),posEci(3));
            thisDate = ephemeris.Epoch;
            thisECI2ECEF = obj.eciFrame.getTransformTo(obj.ecefFrame, thisDate);
            thisPositionECEF = thisECI2ECEF.getRotation().applyTo(thisPositionECI);
            positionECEF = [thisPositionECEF.getX(),thisPositionECEF.getY(),thisPositionECEF.getZ()];
        end
        
        function displayLinkedSatellites(obj,pointLla, starConstellation, ...
                        posIntervals, posSpan, mask, step)
            linkedSatelliteStates = obj.getLinkedSatellites(pointLla, ...
                starConstellation, posIntervals, posSpan, mask);
            NoOfIntervals = size(linkedSatelliteStates,2);
            allStates = linkedSatelliteStates{1,step};
            allEcef = zeros(size(allStates, 1), 3);
            allEcef = allStates(:,2:4);
            if step > NoOfIntervals
                step = NoOfIntervals;
            elseif step < 1 
                step = 1;
            end
            pointEcef = obj.llaToEcef(pointLla);
            geoMap = imread("blue.PNG");
            figure()
            [X,Y,Z] = sphere(50);
            hold on
            rotate3d
            axis equal
            mesh(X.*6378137,Y.*6378137,Z.*6378137,flipud(geoMap),'FaceColor','texturemap',EdgeColor='none')
            plot3(pointEcef(1), pointEcef(2), pointEcef(3), 'r.','MarkerSize',25)
            for i = 1:size(allEcef, 1)
                plot3(allEcef(i, 1), allEcef(i, 2), allEcef(i, 3), 'k.','MarkerSize',25)
            end
            axis off
        end
    
        function lgnssEphemeris = getNavigationMessageFromEstimatorPropagator(obj, ...
                thisPropagator, messageFitInterval)
            
            initialDate = thisPropagator.getMinDate();
            observationTimeStep = 1.0;
            import org.orekit.utils.*
            R_e = Constants.WGS84_EARTH_EQUATORIAL_RADIUS;
            mu = Constants.WGS84_EARTH_MU;
            trueEphemeris = zeros(floor((messageFitInterval/observationTimeStep)) + 1,7);
            trueEphemerisEcef = zeros(floor((messageFitInterval/observationTimeStep)) + 1,7);
            unixTimes = zeros(floor((messageFitInterval/observationTimeStep)) + 1,1);
            currentDateTime = initialDate;
            mn = 0;
            while currentDateTime.compareTo(initialDate.shiftedBy(messageFitInterval)) <=0
                mn = mn + 1;
                trueEphemeris(mn , 1) = mn - 1;
                trueEphemerisEcef(mn , 1) = mn - 1;
                unixTimes(mn , 1) = obj.getUnixTime(currentDateTime);
                truePV = thisPropagator.propagate(currentDateTime).getPVCoordinates();
                trueEphemeris(mn , 2 : 7) = [truePV.getPosition().getX(), truePV.getPosition().getY(), ...
                    truePV.getPosition().getZ(), truePV.getVelocity().getX(), ...
                    truePV.getVelocity().getY(), truePV.getVelocity().getZ()];

                truePVEcef = thisPropagator.propagate(currentDateTime).getPVCoordinates(obj.ecefFrame);
                trueEphemerisEcef(mn , 2 : 7) = [truePVEcef.getPosition().getX(), truePVEcef.getPosition().getY(), ...
                    truePVEcef.getPosition().getZ(), truePVEcef.getVelocity().getX(), ...
                    truePVEcef.getVelocity().getY(), truePVEcef.getVelocity().getZ()];
                currentDateTime = currentDateTime.shiftedBy(observationTimeStep);
            end
            if (mn) < size(trueEphemeris , 1)
                ...
            end
            orbit_data.pos_m = trueEphemeris(:,2:4); %ECI
            orbit_data.vel_m_s = trueEphemeris(:,5:7); %ECI
            orbit_data.elapsed_time_sec = trueEphemeris(:,1);
            ICUIS = obj.ECI2COE(trueEphemeris(1,2:4),trueEphemeris(1,5:7));
            SMA = ICUIS.a; % [m]

            % Fit interval of interest.
            fit_interval = messageFitInterval; % [seconds]

            % File names.
            file_altitudes = SMA - R_e; % [km]

            % Total time in file.
            time_in_file = orbit_data.elapsed_time_sec(end); % [sec]

            % Start the clock.
            tic;

            % Get the analytical coefficients for the URE equations.
            [coeff_A2, coeff_C2, coeff_R2, ~, ~, ~] = ...
                obj.analytic_URE_eqn(  file_altitudes );

            % Determine the message start times.
            message_start_times = ...
                0:messageFitInterval:(time_in_file - fit_interval);

            % Determine the number of messages that can be made per file.
            num_eph_per_file = length( message_start_times );

            % Initialize variables.

            lgnssEphemeris(num_eph_per_file).Asqrt  = [];
            lgnssEphemeris(num_eph_per_file).e      =[];
            lgnssEphemeris(num_eph_per_file).i0     =[];
            lgnssEphemeris(num_eph_per_file).Omega0 =[];
            lgnssEphemeris(num_eph_per_file).Omega  =[];
            lgnssEphemeris(num_eph_per_file).M0     =[];

            lgnssEphemeris(num_eph_per_file).Cus    =[];
            lgnssEphemeris(num_eph_per_file).Cuc    =[];
            lgnssEphemeris(num_eph_per_file).Crs    =[];
            lgnssEphemeris(num_eph_per_file).Crc    =[];
            lgnssEphemeris(num_eph_per_file).Cis    =[];
            lgnssEphemeris(num_eph_per_file).Cic    =[];

            lgnssEphemeris(num_eph_per_file).IDOT      =[];
            lgnssEphemeris(num_eph_per_file).Omega_dot =[];
            lgnssEphemeris(num_eph_per_file).Delta_n   =[];
            lgnssEphemeris(num_eph_per_file).Toe   =[];
            lgnssEphemeris(num_eph_per_file).Epoch = [];

            idx_message = 1;

            % Define the time vector for fitting.
            time = ...
                orbit_data.elapsed_time_sec(...
                orbit_data.elapsed_time_sec <= fit_interval );

            % Get the start / end index for the data to fit to.
            idx_start = find(orbit_data.elapsed_time_sec == ...
                message_start_times(idx_message));
            idx_end = find(orbit_data.elapsed_time_sec == ...
                message_start_times(idx_message) + fit_interval);

            % Get position and velocity vectors
            pos = ...
                orbit_data.pos_m(idx_start:idx_end,:);
            vel = ...
                orbit_data.vel_m_s(idx_start:idx_end,:);

            theta_g = obj.utc2gmst( datevec(unixTimes(1)) ); % [rad]

            % Form initial guess with the 6 Keplerian elements.
            [coe, ~, ~] = obj.ECI2COE(trueEphemeris(1,2:4),trueEphemeris(1,5:7));

            % Form the initial guess for the estimator.
            a = coe.a; % [m]
            n = sqrt(mu/a^3); % [rad/sec]
            ecc = coe.e; % [-]
            inc = coe.i * pi / 180; % [rad]
            RAAN = coe.RAAN * pi / 180; % [rad]
            omega = coe.omega * pi /180; % [rad]
            M0 = coe.M * pi / 180; % [rad]

            Cus = 0; % [rad]
            Cuc = 0; % [rad]
            Crs = 0; % [rad]
            Crc = 0; % [rad]
            Cis = 0; % [rad]
            Cic = 0; % [rad]

            IDOT = 0; % [rad/s]
            OMEGA_DOT = 0; % [rad/s]
            delta_n = 0; % [rad/s]

            % Form the initial guess.
            initial_guess(1) = a;
            initial_guess(2) = ecc;
            initial_guess(3) = inc;
            initial_guess(4) = RAAN;
            initial_guess(5) = omega;
            initial_guess(6) = M0;

            initial_guess(7)  = Cus;
            initial_guess(8)  = Cuc;
            initial_guess(9)  = Crs;
            initial_guess(10) = Crc;
            initial_guess(11) = Cis;
            initial_guess(12) = Cic;

            initial_guess(13) = IDOT;
            initial_guess(14) = OMEGA_DOT;
            initial_guess(15) = delta_n;

            % Define the convergence criteria.
            ConvCrit = 1e-11;

            % Fit ephemeris parameters or subset.
            % Define the weighting matrix, use the identity matrix for now.
            Wmat = ones( size(time) );

            fit_parameters = ones(1,15);

            % Fit parameters.
            [a, ecc, inc, RAAN, omega, M0,...
                Cus, Cuc, Crs, Crc, Cis, Cic, ...
                IDOT, OMEGA_DOT, delta_n, ~, ~, ~] = ...
                obj.COE15_estimator_wrapper(time, pos, vel, initial_guess, ...
                Wmat, ConvCrit, fit_parameters, ...
                theta_g, coeff_R2, coeff_A2, coeff_C2);

            % Error analysis.
            [~, ~, ~, ...
                ~, ~, ~, ...
                lgnssEphemeris] = obj.eph_error_analysis(sqrt(a), ecc, inc, RAAN, omega, M0, ...
                Cus, Cuc, Crc, Crs, Cic, Cis, ...
                IDOT, OMEGA_DOT, delta_n, time, pos, vel, theta_g, ...
                coeff_R2, coeff_A2, coeff_C2);
            lgnssEphemeris(num_eph_per_file).Epoch = initialDate;
        end

        function saveRICError(obj, mainScenarioData, orb, closeStat)
            ricposdiff = obj.getRICPositionDifference(mainScenarioData(orb).realPropagator,...
                mainScenarioData(orb).estimatedPropagator);
            figureInfo = obj.getFigureTitle(mainScenarioData(orb).linkList(1));
            thisFigureName = ['thisFigureCaseNo' num2str(orb) '.fig'];
            thisFigure = obj.getErrorFigure(ricposdiff, figureInfo);
            saveas(thisFigure, thisFigureName)
            if closeStat
                close(thisFigure)
            end
        end

        function positionErrors = getUserPositioningError(obj, ...
                inputStructure, linkedSatelliteStates, mainScenarioData)
            pointLla = inputStructure.underPositioningPoint;
            posSpan = inputStructure.positioningTimeSpanInHours;
            posIntervals = inputStructure.positioningIntervals;
            pointEcef = obj.llaToEcef(pointLla);
            if posSpan*3600.0 > 30.0 * 60.0
                positionErrors = [];
                return
            end
            positionErrors = zeros(posIntervals + 1, 4);

            PsuedoRangeAndEcef = cell(1, posIntervals + 1);
            intervalsDuration = (posSpan*3600.0/(posIntervals + 1));
            for i = 1:posIntervals + 1
                linkedStates = linkedSatelliteStates{1,i};
                if isempty(linkedStates)
                else
                    indexArray = linkedStates(:,1);
                    states = zeros(size(indexArray,1), 4);
                    for j = 1:size(indexArray,1)
                        for z = 1:size(mainScenarioData,2)
                            if indexArray(j) == mainScenarioData(z).underODSatellite
                                EphemerisMessage = mainScenarioData(z).lgnssEphemeris;
                                break
                            end
                        end
                        states(j,1:3) = obj.getEcefFromNavigationMessage(EphemerisMessage, intervalsDuration*(i-1));
                        states(j,4) = norm(states(j,1:3) - pointEcef) - 0.87 + (1.7334)*rand(1,1);
                    end
                    PsuedoRangeAndEcef{1,i} = states;
                end
            end

            for i = 1:posIntervals + 1
                linkedStates = PsuedoRangeAndEcef{1,i};
                if isempty(linkedStates)
                elseif size(linkedStates, 1) < 4
                else
                    positions = linkedStates(:,1:3)'; % u were here
                    ranges = linkedStates(:, 4)';
                    weights = diag(ones(1,length(ranges)));
                    [N1,~] = obj.trilateration(positions, ranges, weights);
                    estimatedPosition = N1(2:4,1)';
                    positionErrors(i,1:3) = estimatedPosition - pointEcef;
                    positionErrors(i,4) = norm(estimatedPosition - pointEcef);
                end
            end


        end

        function  obj = addGroundStation(obj, latitude, longitude, altitude, name)
            import org.orekit.bodies.*;
            import org.orekit.frames.*;
            import org.orekit.estimation.measurements.*;
            stationPoint = GeodeticPoint(latitude * obj.DEG2RAD, longitude * obj.DEG2RAD, altitude);
            stationFrame = TopocentricFrame(obj.wgs84Ellipsoid, stationPoint, name);
            groundStation = GroundStation(stationFrame);
            obj.groundStationsList = [obj.groundStationsList; groundStation];
            obj.numberOfGroundStation = obj.numberOfGroundStation + 1;
        end

        function obj = feedGroundStations(obj, names)
            addedGroundStations = [];
            for i = 1:size(names,2)
                if strcmpi(names(i),"Boushehr")
                    obj = obj.addGroundStation(30.0, 48.5, 0.0, "Boushehr");
                    addedGroundStations = [addedGroundStations; names(i)];
                elseif strcmpi(names(i),"Tehran")
                    obj = obj.addGroundStation(35.0, 51.0, 0.0, "Tehran");
                    addedGroundStations = [addedGroundStations; names(i)];
                elseif strcmpi(names(i),"Aleshtar")
                    obj = obj.addGroundStation(33.0, 48.0, 1500.0, "Aleshtar");
                    addedGroundStations = [addedGroundStations; names(i)];
                elseif strcmpi(names(i),"Bandar Abbas")
                    obj = obj.addGroundStation(25.0, 59.0, 0.0, "Bandar Abbas");
                    addedGroundStations = [addedGroundStations; names(i)];
                elseif strcmpi(names(i),"Mashhad")
                    obj = obj.addGroundStation(35.0, 58.0, 0.0, "Mashhad");
                    addedGroundStations = [addedGroundStations; names(i)];
                elseif strcmpi(names(i),"Tabriz")
                    obj = obj.addGroundStation(38.0, 46.0, 0.0, "Tabriz");
                    addedGroundStations = [addedGroundStations; names(i)];
                elseif strcmpi(names(i),"Sistan")
                    obj = obj.addGroundStation(27.0, 63.0, 0.0, "Sistan");
                    addedGroundStations = [addedGroundStations; names(i)];
                else
                    pass
                end
                % print(" ")
                % print("added Ground Stations ({}): {}".format(len(addedGroundStations), addedGroundStations))

            end
        end

    end

end