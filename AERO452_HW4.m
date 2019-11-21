%% AERO452-01 HW 6
%% Lacey Davis 
% Due Date: December 3, 2019

clear all; close all; clc;

%constants:
mu_e = 398600 ; %km3/s2
Re = 6378 ; %km
Rse = 149600000 ; %km 
SF = 1367 ; %W/m2
Psr = 4.57*10^-6 ; %N/m2
c = 3*10^8 ; %m/s

%% Question 1
disp('Question 1: ')

R1 = [-26175.1034, 12757.0706, 14626.6556] ; 
V1 = [2.376641, .139677, 2.078097] ; 
state1 = [R1 V1] ; 
tspan = 13*24*3600 ; %sec
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;
[dstate1] = ode45(@linearmotion, [0 tspan], state1, options) ;

R = dstate1.y(1:3,:) ; 

for ii = 1:length(dstate1.x) 
    
    jj = 10:((23-10)/length(dstate1.x)):23 ; 
    JD(ii) = julian(8, jj(ii), 1997, 0, 0) ;
    [rsun1(1:3,ii),rtasc1(ii),decl1(ii)] = sun ( JD(ii) ) ;
    rsun1(1:3,ii) = rsun1(1:3,ii)*Rse ; %AU to km
    
    r1(ii) = norm(R(1:3,ii)) ; 
    rs(ii) = norm(rsun1(1:3,ii)) ;
    theta(ii) = acos((dot(rsun1(1:3,ii),R(1:3,ii)))/(rs(ii)*r1(ii))) ;  
    theta1(ii) = acos(Re/(r1(ii))) ; 
    theta2(ii) = acos(Re/(rs(ii))) ;
    xx(ii) = (theta1(ii) + theta2(ii)) ; 
        if xx(ii) <= theta(ii) 
            F(ii) = 0 ; 
        else 
            F(ii) = 1 ; 
        end 
end

off = find(~F) ; 
ecl = 10+(off*((23-10)/length(dstate1.x))) ;
%     for mm = 1:length(off) 
%         disp('Day.Time Fraction of Eclipse from Aug 10 to Aug 23: ')
%         disp(jj(off(mm)))
%     end 
per_ecl = (length(dstate1.x)-length(ecl))/length(dstate1.x) ; 

figure(1) 
plot((dstate1.x/(24*3600)),F)

disp('Time in Eclipse is shown by the vertical spikes in the figure. The cycle of eclipse to sun exposure is cyclical as expected for each day.')
disp('The spacecraft spends ')
disp(per_ecl*100) 
disp('percent of the time span exposed to the Sun where SRP has perturbational effects. ')
disp('If the eclipse time is ignored -- the absence of the effect of the shadow function (which is a linear relationship to SRP) ')
disp('would increase the SRP by ')
disp(100 - (per_ecl*100))
disp('percent. This is significant enough to not ignore Earth eclipse when performing calculations. ')

%% Question 2
disp('Question 2: ') 
disp('Initial Julian Date from Exmaple 12.11 in Curtis. ')

%HEO Orbit Graphs 
tspan2 = 60*24*3600 ; %sec
h2 = 69084.1 ; 
ecc2 = .741 ;
RAAN2 = deg2rad(0) ; 
inc2 = deg2rad(63.4) ;
arg2 = deg2rad(270) ; 
theta2 = deg2rad(0) ; 
a2 = 26553.4 ; 
T2 = 11.9616*3600 ; %sec
state2 = [h2, inc2, RAAN2, arg2, ecc2, theta2] ; 
[tcoes dcoes] = ode45(@vop, [0 tspan2], state2, options) ; 

    for nn = 1:length(tcoes) 
        a2n = dcoes(nn,1)^2/mu_e/(1 - dcoes(nn,5)^2);
        ra2n(nn) = a2n*(1+dcoes(nn,5)) ; 
        rp2n(nn) = a2n*(1-dcoes(nn,5)) ;
    end 

figure(2) 
subplot(3,1,1)
plot((tcoes/(24*3600)),rad2deg(dcoes(:,3)))
title('HEO: Change in RAAN')
xlabel('Time (Days)')
ylabel('Degrees') 
subplot(3,1,2)
plot((tcoes/(24*3600)),rad2deg(dcoes(:,2)))
title('HEO: Change in Inclination')
xlabel('Time (Days)')
ylabel('Degrees')
subplot(3,1,3)
plot((tcoes/(24*3600)),rad2deg(dcoes(:,4)))
title('HEO: Change in Argument of Perigee')
xlabel('Time (Days)')
ylabel('Degrees')

%% Question 3
disp('Question 3: ') 

%Given: 
aa = [-1.05, 0] ;
av = [0, -0.066429] ;
at = 2*pi ; 

bb = [-1.15, 0] ;
bv = [0, -0.0086882909] ;
bt = 3.6 ;

cc = [.1, 0] ;
cv = [-3.35, 3] ;
ct = 3.6 ;

dd = [.1, 0] ; 
dv = [-3.37, 3] ;
dt = 6 ;

ee = [.1, 0] ;
ev = [-3.4, 3] ;
et = 6 ;

ff = [1.25, 0] ;
fv = [0, 0] ;
ft = 2*pi ;

gg = [-.5, 0] ;
gv = [0, 0] ;
gt = 2*pi ;

hh = [-1.1, 0] ;
hv = [0, 0] ;
ht = 2*pi ;



% Functions 

    %Julian Date
      function [ JD ] = julian( m,d,y,hr,min )
%julian is a function that will convert UT into Julian Date
%   inputs are month, day, year, and UT time fraction

Jo = 367*y - floor((7*(y+floor((m+9)/12)))/4) + floor((275*m)/9) + d + 1721013.5 ;
tf = hr+(min/60) ; 
JD = Jo + (tf/24) ;

      end

      %Two Body Motion 
function dstatedt = linearmotion (t, state)
%function for ode45 proces, defines the differential functions to integrate
mu_e = 398600 ;
R = [state(1) state(2) state(3)] ; 
V = [state(4) state(5) state(6)] ;
r = norm([state(1) state(2) state(3)]) ; %norm of the position vector
h = norm(cross(R,V)) ;

dx = state(4) ; %velocity differential equations
dy = state(5) ;
dz = state(6) ;

r = norm([state(1) state(2) state(3)]) ;    %norm of the position vector

ddx = (-mu_e * state(1)) / r^3 ;  %Equations of relative motion 
ddy = (-mu_e * state(2)) / r^3 ;
ddz = (-mu_e * state(3)) / r^3 ;

dstatedt = [dx;dy;dz;ddx;ddy;ddz] ;    

end

%VoP 
function [dcoes] = vop(t, state) 
% u is some osculating element, use ODE45
 mu = 398600 ; 
 mu_s = 132712440018 ; 
 Re = 6378 ; 
    h = state(1) ; 
    inc = state(2) ; 
    RAAN = state(3) ; 
    arg = state(4) ;
    ecc = state(5) ; 
    theta = state(6) ; 
    
    coes = [state(1) state(5) state(3) state(2) state(4) state(6)] ;
 %Coes to Rv
 [R, V] = sv_coes(coes,mu) ; 
 r = norm(R) ; 
 H = cross(R,V) ;
 
 alt = r-Re ;


function [lookfor stop direction] = terminate(t,y)

% From Curtis:

%

% This function specifies the event at which ode45 terminates.

% ????????????????????????????????????????????????????????????????????????
lookfor = alt - 100; % = 0 when altitude = 100 km

stop = 1; % 1 means terminate at lookfor = 0; Otherwise 0

direction = -1; % -1 means zero crossing is from above

 end   
 
    uR = R/norm(R) ; 
    uN = H/norm(H) ;
    uT = cross(uN, uR) ; 
    
%Julian Date Calculation   
JDo = 2454283 ; %initial start time from example 12.12
tf = t/(24*3600) ; 
JDf = JDo + tf ; %JDf will be 60 days from JDo

%Sun Position 
Rse = [0;149600000;0] ; %km 
[rsun,~,~] = sun (JDf) ;
rsun = rsun'*norm(Rse) ; %AU to km

Rss = rsun - R ;
q = dot(R,((2*rsun)-R))/(norm(rsun)^2) ; 
F = (((q^2) - (3*q) + 3)/(1 + (1 - q)^(3/2))) * q ; 
ap = (mu_s*((F*rsun)-R))/(norm(Rss)^3) ;
    
    %Getting ap in RTN vector form 
    aR = dot(ap,uR) ; 
    aN = dot(ap,uN) ;
    aT = dot(ap,uT) ; 
    
    %magnitudes of COEs
    dh = r*aT ; 
    
    decc = ((h/mu)*aR*sin(theta)) + ...
        (((1/(mu*h))*(((h^2)+(mu*r)*cos(theta)) + (mu*ecc*r)))*aT) ;
    
    dtheta = (h/r^2) + ...
        ((1/(ecc*h))*((((h^2)/mu)*aR*cos(theta)) - ...
        ((((h^2)/mu)+r)*aT*sin(theta)))) ;
    
    dinc = (r/norm(h))*aN*cos(arg+theta) ; 
    
    dRAAN = ((r*sin(arg+theta))/(norm(h)*sin(inc))) * aN ; 
    
    darg = (((-r*sin(arg+theta))/(h*tan(inc)))*aN) - ...
        ((1/(ecc*h))*((((h^2)/mu)*aR*cos(theta)) - ...
        ((((h^2)/mu)+r)*aT*sin(theta)))) ; 
    
  dcoes = [dh;dinc;dRAAN;darg;decc;dtheta] ;  
    
end  

%Coes to RV
function [r, v] = sv_coes(coe,mu)
h = coe(1);

e = coe(2);

RA = coe(3);

incl = coe(4);

w = coe(5);

TA = coe(6);

%...Equations 4.45 and 4.46 (rp and vp are column vectors):

rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);

vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

%...Equation 4.34:

R3_W = [ cos(RA) sin(RA) 0

 -sin(RA) cos(RA) 0

 0 0 1];

%...Equation 4.32:

R1_i = [1 0 0

 0 cos(incl) sin(incl)

 0 -sin(incl) cos(incl)];

%...Equation 4.34:

R3_w = [ cos(w) sin(w) 0

 -sin(w) cos(w) 0

 0 0 1];

%...Equation 4.49:

Q_pX = (R3_w*R1_i*R3_W)';

%...Equations 4.51 (r and v are column vectors):

r = Q_pX*rp;

v = Q_pX*vp;

end

      %Planet Position -- adapted from Dr. A
function [planet_coes] = AERO451planetary_elements2(planet_id,T)
% Planetary Ephemerides from Meeus (1991:202-204) and J2000.0
% Output:
% planet_coes
% a = semimajor axis (km)
% ecc = eccentricity
% inc = inclination (degrees)
% raan = right ascension of the ascending node (degrees)
% w_hat = longitude of perihelion (degrees)
% L = mean longitude (degrees)

% Inputs:
% planet_id - planet identifier:
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 6 = Saturn
% 7 = Uranus
% 8 = Neptune

if planet_id == 1
    a = 0.387098310; % AU but in km later
    ecc = 0.20563175 + 0.000020406*T - 0.0000000284*T^2 - 0.00000000017*T^3;
    inc = 7.004986 - 0.0059516*T + 0.00000081*T^2 + 0.000000041*T^3; %degs
    raan = 48.330893 - 0.1254229*T-0.00008833*T^2 - 0.000000196*T^3; %degs
    w_hat = 77.456119 +0.1588643*T -0.00001343*T^2+0.000000039*T^3; %degs
    L = 252.250906+149472.6746358*T-0.00000535*T^2+0.000000002*T^3; %degs
elseif planet_id == 2
    a = 0.723329820; % AU
    ecc = 0.00677188 - 0.000047766*T + 0.000000097*T^2 + 0.00000000044*T^3;
    inc = 3.394662 - 0.0008568*T - 0.00003244*T^2 + 0.000000010*T^3; %degs
    raan = 76.679920 - 0.2780080*T-0.00014256*T^2 - 0.000000198*T^3; %degs
    w_hat = 131.563707 +0.0048646*T -0.00138232*T^2-0.000005332*T^3; %degs
    L = 181.979801+58517.8156760*T+0.00000165*T^2-0.000000002*T^3; %degs
elseif planet_id == 3 
    a = 1.000001018; % AU
    ecc = 0.01670862 - 0.000042037*T - 0.0000001236*T^2 + 0.00000000004*T^3;
    inc = 0.0000000 + 0.0130546*T - 0.00000931*T^2 - 0.000000034*T^3; %degs
    raan = 0.0; %degs
    w_hat = 102.937348 + 0.3225557*T + 0.00015026*T^2 + 0.000000478*T^3; %degs
    L = 100.466449 + 35999.372851*T - 0.00000568*T^2 + 0.000000000*T^3; %degs
elseif planet_id == 4
    a = 1.523679342; % AU
    ecc = 0.09340062 + 0.000090483*T - 0.00000000806*T^2 - 0.00000000035*T^3;
    inc = 1.849726 - 0.0081479*T - 0.00002255*T^2 - 0.000000027*T^3; %degs
    raan = 49.558093 - 0.2949846*T-0.00063993*T^2 - 0.000002143*T^3; %degs
    w_hat = 336.060234 +0.4438898*T -0.00017321*T^2+0.000000300*T^3; %degs
    L = 355.433275+19140.2993313*T+0.00000261*T^2-0.000000003*T^3; %degs
elseif planet_id == 5
    a = 5.202603191 + 0.0000001913*T; % AU
    ecc = 0.04849485+0.000163244*T - 0.0000004719*T^2 + 0.00000000197*T^3;
    inc = 1.303270 - 0.0019872*T + 0.00003318*T^2 + 0.000000092*T^3; %degs
    raan = 100.464441 + 0.1766828*T+0.00090387*T^2 - 0.000007032*T^3; %degs
    w_hat = 14.331309 +0.2155525*T +0.00072252*T^2-0.000004590*T^3; %degs
    L = 34.351484+3034.9056746*T-0.00008501*T^2+0.000000004*T^3; %degs
elseif planet_id == 6
    a = 9.5549009596 - 0.0000021389*T; % AU
    ecc = 0.05550862 - 0.000346818*T -0.0000006456*T^2 + 0.00000000338*T^3;
    inc = 2.488878 + 0.0025515*T - 0.00004903*T^2 + 0.000000018*T^3; %degs
    raan = 113.665524 - 0.2566649*T-0.00018345*T^2 + 0.000000357*T^3; %degs
    w_hat = 93.056787 +0.5665496*T +0.00052809*T^2-0.000004882*T^3; %degs
    L = 50.077471+1222.1137943*T+0.00021004*T^2-0.000000019*T^3; %degs
elseif planet_id == 7
    a = 19.218446062-0.0000000372*T+0.00000000098*T^2; % AU
    ecc = 0.04629590 - 0.000027337*T + 0.0000000790*T^2 + 0.00000000025*T^3;
    inc = 0.773196 - 0.0016869*T + 0.00000349*T^2 + 0.00000000016*T^3; %degs
    raan = 74.005947 + 0.0741461*T+0.00040540*T^2 +0.000000104*T^3; %degs
    w_hat = 173.005159 +0.0893206*T -0.00009470*T^2+0.000000413*T^3; %degs
    L = 314.055005+428.4669983*T-0.00000486*T^2-0.000000006*T^3; %degs
elseif planet_id == 8
    a = 30.110386869-0.0000001663*T+0.00000000069*T^2; % AU
    ecc = 0.00898809 + 0.000006408*T -0.0000000008*T^2;
    inc = 1.769952 +0.0002557*T +0.00000023*T^2 -0.0000000000*T^3; %degs
    raan = 131.784057 - 0.0061651*T-0.00000219*T^2 - 0.000000078*T^3; %degs
    w_hat = 48.123691 +0.0291587*T +0.00007051*T^2-0.000000000*T^3; %degs
    L = 304.348665+218.4862002*T+0.00000059*T^2-0.000000002*T^3; %degs
end

planet_coes = [a;ecc;inc;raan;w_hat;L];
%Convert to km:
au = 149597870;
planet_coes(1) = planet_coes(1)*au;
end
      %Sun Position -- adapted from Vallado
%
% ------------------------------------------------------------------------------
%
%                           function sun
%
%  this function calculates the geocentric equatorial position vector
%    the sun given the julian date.  this is the low precision formula and
%    is valid for years from 1950 to 2050.  accuaracy of apparent coordinates
%    is 0.01  degrees.  notice many of the calculations are performed in
%    degrees, and are not changed until later.  this is due to the fact that
%    the almanac uses degrees exclusively in their formulations.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%    vallado     - fix mean lon of sun                            7 mat 2004
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    rsun        - ijk position vector of the sun au
%    rtasc       - right ascension                rad
%    decl        - declination                    rad
%
%  locals        :
%    meanlong    - mean longitude
%    meananomaly - mean anomaly
%    eclplong    - ecliptic longitude
%    obliquity   - mean obliquity of the ecliptic
%    tut1        - julian centuries of ut1 from
%                  jan 1, 2000 12h
%    ttdb        - julian centuries of tdb from
%                  jan 1, 2000 12h
%    hr          - hours                          0 .. 24              10
%    min         - minutes                        0 .. 59              15
%    sec         - seconds                        0.0  .. 59.99          30.00
%    temp        - temporary variable
%    deg         - degrees
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 281, alg 29, ex 5-1
%
% [rsun,rtasc,decl] = sun ( jd );
% ------------------------------------------------------------------------------

function [rsun,rtasc,decl] = sun ( jd )

        twopi      =     2.0*pi;
        deg2rad    =     pi/180.0;

        % -------------------------  implementation   -----------------
        % -------------------  initialize values   --------------------
        tut1= ( jd - 2451545.0  )/ 36525.0;
%fprintf(1,'tut1 %14.9f \n',tut1);

        meanlong= 280.460  + 36000.77*tut1;
        meanlong= rem( meanlong,360.0  );  %deg

        ttdb= tut1;
        meananomaly= 357.5277233  + 35999.05034 *ttdb;
        meananomaly= rem( meananomaly*deg2rad,twopi );  %rad
        if ( meananomaly < 0.0  )
            meananomaly= twopi + meananomaly;
        end

        eclplong= meanlong + 1.914666471 *sin(meananomaly) ...
                    + 0.019994643 *sin(2.0 *meananomaly); %deg
        eclplong= rem( eclplong,360.0  );  %deg

        obliquity= 23.439291  - 0.0130042 *ttdb;  %deg

        eclplong = eclplong *deg2rad;
        obliquity= obliquity *deg2rad;

        % --------- find magnitude of sun vector, )   components ------
        magr= 1.000140612  - 0.016708617 *cos( meananomaly ) ...
                              - 0.000139589 *cos( 2.0 *meananomaly );    % in au's

        rsun(1)= magr*cos( eclplong );
        rsun(2)= magr*cos(obliquity)*sin(eclplong);
        rsun(3)= magr*sin(obliquity)*sin(eclplong);

%fprintf(1,'meanlon %11.6f meanan %11.6f eclplon %11.6f obli %11.6f \n', ...
%           meanlong,meananomaly/deg2rad,eclplong/deg2rad,obliquity/deg2rad);
%fprintf(1,'rs %11.9f %11.9f %11.9f \n',rsun);
%fprintf(1,'magr %14.7f \n',magr);

        rtasc= atan( cos(obliquity)*tan(eclplong) );

        % --- check that rtasc is in the same quadrant as eclplong ----
        if ( eclplong < 0.0  )
            eclplong= eclplong + twopi;    % make sure it's in 0 to 2pi range
        end
        if ( abs( eclplong-rtasc ) > pi*0.5  )
            rtasc= rtasc + 0.5 *pi*round( (eclplong-rtasc)/(0.5 *pi));
        end
        decl = asin( sin(obliquity)*sin(eclplong) );
end