%% Jupiter and Venus N-Body 

%inputs: initial Julian Date, initial COEs, initial state, 
%outputs: R, V, coes

%functions needed: state to coes, coes to sv, enckesnb, uni-variable, position of Sun, position of planets

%propagation method: enckes, ensures that time step is converted to days to
%find position of the bodies 

clear all; close all; clc;

%constants:
mu_e = 398600 ; %km3/s2

%Hw4 Test: 
%HEO Orbit Graphs 
% h = 69084.1 ; 
% ecc = .741 ;
% RAAN = deg2rad(2*pi) ; 
% inc = deg2rad(63.4) ;
% arg = deg2rad(270) ; 
% theta = deg2rad(0) ; 
% a = 26553.4 ; 
% T = 11.9616*60 ; %sec

coes = [h ecc RAAN inc arg theta] ;
tf = 60*24*3600 ; 
[r, v] = sv_coes(coes,mu_e) ;
dt = 300 ; %time step in seconds 
[R_enc, V_enc, tt] = enckesnb(r, r, v, v, [0;0;0], [0;0;0], dt, tf, mu_e)  ;
    for mm = 1:length(R_enc(1,:)) 
      [henc(mm), incenc(mm), RAANenc(mm), argenc(mm), eccenc(mm), raenc(mm), rpenc(mm)] = state2coes(R_enc(1:3,mm), V_enc(1:3,mm), mu_e) ;
    end 

% figure(1) 
% subplot(1,3,1)
% plot((tt/(24*3600)),(rad2deg(RAANenc)-rad2deg(RAAN)))
% title('HEO: Change in RAAN')
% xlabel('Time (Days)')
% ylabel('Degrees')
% axis([0 60 -.25 0])
% subplot(1,3,2)
% plot((tt/(24*3600)),rad2deg(incenc)-rad2deg(inc))
% title('HEO: Change in Inclination')
% xlabel('Time (Days)')
% ylabel('Degrees')
% axis([0 60 -.005 .025])
% subplot(1,3,3)
% plot((tt/(24*3600)),rad2deg(argenc)-rad2deg(arg))
% title('HEO: Change in Argument of Perigee')
% xlabel('Time (Days)')
% ylabel('Degrees')
% axis([0 60 -.02 .12])

    %Enckes Method With Perturbations -- N Body
function [R, V, tt] = enckesnb(R, R_osc, V, V_osc, dro, dvo, dt, tf, mu_e)  
ii = 1; 
tt(ii) = 0 ;  
rectify = 1 ; 
while ii < 10E5 && tt(ii) < tf
    
        R_osc(:,ii) = R(:,ii) ;
        V_osc(:,ii) = V(:,ii) ; 
        dr(:,ii) = [0; 0; 0] ;
        dv(:,ii) = [0; 0; 0] ; 
         
    r(ii) = norm(R(:,ii)) ; 
    v(ii) = norm(V(:,ii)) ; 
    r_osc(ii) = norm(R_osc(:,ii)) ; 
    v_osc(ii) = norm(V_osc(:,ii)) ; 
    dr(:,1) = dro ; 
    dv(:,1) = dvo ; 

        %Julian Date Calculation   
    JDo = 2454283 ; %initial start time from example 12.12
    JDf = JDo + (tt(ii)/(24*3600)) ; %JDf will be 60 days from JDo

    %Sun Position 
    mu_s = 132712440018 ; 
    Re = 6378 ; 
    Rse = [0;149600000;0] ; %km 
    [rsun,~,~] = sun (JDf) ;
    rsun = rsun'*norm(Rse) ; %AU to km
    
    %Jupiter Position 
    mu_J = 126686534 ; 
    RJs = [0;778300000;0] ; 
    RJe = RJs - Rse ; 
    [jcoes] = AERO451planetary_elements2(5,JDf) ;
    rJ = [0;jcoes(1);0] ;
    
    %Venus Position
    mu_V = 324859 ;
    RVs = [0;108200000;0] ;
    RVe = RVs - Rse ;
    [vcoes] = AERO451planetary_elements2(2,JDf) ;
    rV = [0;vcoes(1);0] ;
    
    Rss = rsun - R(:,ii) ;
    RssJ = rJ - R(:,ii) ;
    RssV = rV - R(:,ii) ;
    
    qs = dot(R(:,ii),((2*rsun)-R(:,ii)))/(norm(rsun)^2) ; 
    Fs = (((qs^2) - (3*qs) + 3)/(1 + (1 - qs)^(3/2))) * qs ; 
    
    %Jupiter
    qJ = dot(R(:,ii),((2*rJ)-R(:,ii)))/(norm(rJ)^2) ; 
    FJ = (((qJ^2) - (3*qJ) + 3)/(1 + (1 - qJ)^(3/2))) * qJ ;
    
    %Venus 
    qV = dot(R(:,ii),((2*rV)-R(:,ii)))/(norm(rV)^2) ; 
    FV = (((qV^2) - (3*qV) + 3)/(1 + (1 - qV)^(3/2))) * qV ;
    
    ap(:,ii) = (mu_s*((Fs*rsun)-R(:,ii)))/(norm(Rss)^3) + ...
        (mu_J*((FJ*rJ)-R(:,ii)))/(norm(RssJ)^3) + ...
        (mu_V*((FV*rV)-R(:,ii)))/(norm(RssV)^3) ;
 
    q = dot(dr(:,ii),((2*R(:,ii))-dr(:,ii)))/(r(ii)^2) ; 
    F = (((q^2) - (3*q) + 3)/(1 + (1 - q)^(3/2))) * q ; 

    da(:,ii) = ((-mu_e*(dr(:,ii) - (F*R(:,ii))))/(r_osc(ii)^3)) + ap(:,ii) ;
    dv(:,ii+1) = (da(:,ii)*dt) + dv(:,ii) ; 
    dr(:,ii+1) = (.5*da(:,ii)*dt^2) + (dv(:,ii)*dt) + dr(:,ii) ;
 
    [R_osc(1:3,ii+1),V_osc(1:3,ii+1)] = univariable(R_osc(1:3,ii),V(1:3,ii),dt,mu_e)   ;

 tt(ii+1) = tt(ii) + dt ;
 R(:,ii+1) = R_osc(:,ii+1) + dr(:,ii+1) ; 
 V(:,ii+1) = V_osc(:,ii+1) + dv(:,ii+1) ; 
 
 alt = norm(R(:,ii+1))-Re ; 
 
 if alt < 100
     disp('Re-entry Due to Drag After ')
     disp(tt(ii)/(24*3600))
     break  
 end 

ii = ii + 1 ; 
 
end
end

%Universal Variable
function [r1,v1] = univariable(r0,v0,dt,mu)

ii = 0; % initiate counter
ratio = 1; % starting out ratio to just have one
alpha = 2/norm(r0) - (norm(v0)^2)/mu; % inverse of semimajor axis
x = sqrt(mu)*abs(alpha)*dt; % initial guess of x
r = norm(r0); % radius magnitude
vr = dot(v0,r0)/r; % magnitude of radial velocity
while abs(ratio)> 10^(-8) && ii<10000 % tolerance and infinite loop stopper
    ii = ii + 1; % counter increase
    z = alpha*x^2; % current run through's z-value
    % Stumpff Equations
    if z > 0 % ellipses
        S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
        C = (1 - cos(sqrt(z)))/z;
    elseif z < 0 % hyperbolas
        S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
        C = (cosh(sqrt(-z)) - 1)/(-z);
    else % parabolas
        S = 1/6;
        C = 1/2;
    end
    % Curtis for Newton's
    f = r*vr/sqrt(mu)*x^2*C + (1 - alpha*r)*x^3*S + r*x - sqrt(mu)*dt;
    fp = r*vr/sqrt(mu)*x*(1 - alpha*x^2*S) + (1 - alpha*r)*x^2*C + r;
    ratio = f/fp; % Newton's Iteration ratio
    x = x - ratio; % New x value based upon ratio
end
% Solve for new position vector
f = 1 - x^2/r*C; % Lagrange co
g = dt - 1/sqrt(mu)*x^3*S; % Lagrange co
r1 = f*r0 + g*v0; % km, new position vector
% Solve for new velocity vector
fdot = sqrt(mu)/r/norm(r1)*(z*S - 1)*x; % Lagrange co
gdot = 1 - x^2/norm(r1)*C; % Lagrange co
v1 = fdot*r0 + gdot*v0; % km/s, new velocity vector

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

%State to COEs
function [h, inc, RAAN, arg, ecc, ra, rp] = state2coes(dRdt, dVdt, mu_e)
eps = 1.e-10;
drdt = norm(dRdt);
dvdt = norm(dVdt);
dvrdt = dot(dRdt,dVdt)/drdt;
H = cross(dRdt,dVdt);
h = norm(H);

%Inclination
inc = acos(H(3)/h);
N = cross([0 0 1],H);
n = norm(N);

%Right Ascension of Ascending Node
    if n ~= 0
     RAAN = acos(N(1)/n);
         if N(2) < 0
         RAAN = 2*pi - RAAN;
         end
    else
     RAAN = 0;
    end

%Eccentricity
E = 1/mu_e*((dvdt^2 - mu_e/drdt)*dRdt - drdt*dvrdt*dVdt);
ecc = norm(E);

%Arg of Perigee
if n ~= 0
         if ecc > eps
         arg = acos(dot(N,E)/n/ecc);
             if E(3) < 0
             arg = 2*pi - arg;
             end
         else
         arg = 0;
         end
else
 arg = 0;
end

%True Anomaly
    if ecc > eps
       theta = acos(dot(E,dRdt)/ecc/drdt);
         if dvrdt < 0
         theta = 2*pi - theta;
         end
    else

 cp = cross(N,dRdt);
     if cp(3) >= 0
     theta = acos(dot(N,dRdt)/n/drdt);
     else
     theta = 2 * pi - acos(dot(N,dRdt)/n/drdt);
     end
    end

%Semimajor Axis
a = h^2/mu_e/(1 - ecc^2);
ra = a*(1+ecc) ; 
rp = a*(1-ecc) ;
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

planet_coes = [a;ecc;raan;inc;w_hat;L];
%Convert to km:
au = 149597870;
planet_coes(1) = planet_coes(1)*au;
end

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
