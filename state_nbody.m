%% Notes to self

%Add SRP Perterbations

%% Part 1 - Propagation

clear all; close all; clc; 
mu_e = 398600 ; 
W = [0 0 72.9211e-6]';  %angular velocity of the earth (rad/s)
r_e = 6378;             %radius of the Earth (km)
tspan = 100*24*3600 ; %sec
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

%% Iridium - Drag, N-Body, Oblateness, SRP
RI = [-4351.6; -5661.3; 15.3] ;
VI = [.3814; -.2730; 7.4567] ;
stateI = [RI VI] ;
A_i = 19.25;   %m^2
m_i = 860;  %kg

%(1)N-Body (2)Drag (3)J2 (4)J3 (5)J4 (6)J5 (7)J6 (8)SRP
onoff = [1 1 0 1 1 1 1 1];

dt = 100 ; %time step in seconds 
[RI2, VI2, tt, dcoesI] = enckesnb(RI, RI, VI, VI, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_i, m_i, r_e)  ;

figure(1)
earth_sphere
hold on
plot3(RI2(1,:),RI2(2,:),RI2(3,:))
grid on
title('Iridium Propagation')

%% COSMOS - N-body, Oblateness, SRP 
RC = [-10401; 39485; 10622] ;
VC = [-2.9744; -.7741; -.0349] ;
stateC = [RC VC] ;
A_c = 40.48;  %m^2
m_c = 1700;  %kg

%(1)N-Body (2)Drag (3)J2 (4)J3 (5)J4 (6)J5 (7)J6 (8)SRP
onoff = [1 0 1 1 1 1 1 1];

dt = 100 ; %time step in seconds 
[RC2, VC2, tt, dcoesC] = enckesnb(RC, RC, VC, VC, [0;0;0], [0;0;0], dt, ...
    tspan, mu_e, onoff, W, A_c, m_c, r_e)  ;

figure(2)
earth_sphere
hold on
plot3(RC2(1,:),RC2(2,:),RC2(3,:))
grid on
title('Iridium Propagation')

%% RADARSAT-2 - N-Body, Oblateness, SRP, Drag
RR = [-2213.2 88.2 6819.2]';
VR = [-6.0621 3.8439 -2.0172]';

stateR = [RR VR] ;
A_r = 15.18;  %m^2
m_r = 2200;  %kg

%(1)N-Body (2)Drag (3)J2 (4)J3 (5)J4 (6)J5 (7)J6 (8)SRP
onoff = [1 1 0 1 1 1 1 1 1];

dt = 100 ; %time step in seconds 
[RR2, VR2, tt, dcoesR] = enckesnb(RR, RR, VR, VR, [0;0;0], [0;0;0], dt, ...
    tspan, mu_e, onoff, W, A_r, m_r, r_e)  ;

figure(3)
earth_sphere
hold on
plot3(RR2(1,:),RR2(2,:),RR2(3,:))
grid on
title('Iridium Propagation')


%% Functions: 
    %Enckes Method With Perturbations -- N Body -- Drag -- J2-6
function [R, V, tt, dcoes] = enckesnb(R, R_osc, V, V_osc, dro, dvo, dt, tf, ...
    mu_e, onoff, omega_e, A_sc, m, r_e)  
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

        %J Constants
    J2 = 0.00108263;
    J3 = (-2.53266e-6)*J2;
    J4 = (-1.619e-6)*J2;
    J5 = (-2.273e-7)*J2;
    J6 = (5.407e-7)*J2; 
    
        %Julian Date Calculation   
    JDo = 2451838 ; %initial start time 
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
    
        %Drag Perterbations
    vrel = (dr(:,ii) - cross(omega_e,R_osc(:,ii)))*10^3;   %relative velocity (m/s)
    CD = 2.2;                               %coefficient of drag

    altell = (r_osc(ii)-6378);        %elliptical altitude (m)
    rho = atmosphere(altell);       %kg/m^3

    a_drag(:,ii) = -(1/2)*CD*(A_sc/m)*rho*(norm(vrel)^2)*vrel/norm(vrel); %ECI (m/s^2)
    a_drag(:,ii) = a_drag(:,ii)*10^-3;      %km/s^2
    
        %J2-J6 Perterbations
    %J2 perterbations
    a_J2(1,ii) = -(3*J2*mu_e*r_e^2*R_osc(1,ii))/(2*r_osc(ii)^5)*(1 - ...
        ((5*R_osc(3,ii)^2)/r_osc(ii)^2));
    a_J2(2,ii) = -(3*J2*mu_e*r_e^2*R_osc(2,ii))/(2*r_osc(ii)^5)*(1 - ...
        ((5*R_osc(3,ii)^2)/r_osc(ii)^2));
    a_J2(3,ii) = -(3*J2*mu_e*r_e^2*R_osc(3,ii))/(2*r_osc(ii)^5)*(3 - ...
        ((5*R_osc(3,ii)^2)/r_osc(ii)^2));
%     temp = -((3*J2*mu_e*(r_e^2))/(2*r_osc(ii)^5)); 
%     aI = temp*R_osc(1,ii)*(1 - ((5*R_osc(3,ii)^2)/(r_osc(ii)^2)));
%     aJ = temp*R_osc(2,ii)*(1 - ((5*R_osc(3,ii)^2)/(r_osc(ii)^2)));
%     aK = temp*R_osc(3,ii)*(3 - ((5*R_osc(3,ii)^2)/(r_osc(ii)^2)));
%     a_J2(:,ii) = [aI , aJ , aK]';

    %J3 perterbations
    a_J3(1,ii) = -(5*J3*mu_e*r_e^3*R_osc(1,ii)/(2*r_osc(ii)^7))*(3*R_osc(3,ii) ...
        - (7*R_osc(3,ii)^3)/r_osc(ii)^2);
    a_J3(2,ii) = -(5*J3*mu_e*r_e^3*R_osc(2,ii)/(2*r_osc(ii)^7))*(3*R_osc(3,ii) ...
        - (7*R_osc(3,ii)^3)/r_osc(ii)^2);
    a_J3(3,ii) = -(5*J3*mu_e*r_e^3/(2*r_osc(ii)^7))*(6*R_osc(3,ii)^2 ...
        - ((7*R_osc(3,ii)^4)/r_osc(ii)^2) - (3*r_osc(ii)^2)/5);

    %J4 perterbations
    a_J4(1,ii) = (15*J4*mu_e*r_e^4*R_osc(1,ii)/(8*r_osc(ii)^7))*(1-...
        (14*R_osc(3,ii)^2)/r_osc(ii)^2 + (21*R_osc(3,ii)^4)/r_osc(ii)^4);
    a_J4(2,ii) = (15*J4*mu_e*r_e^4*R_osc(2,ii)/(8*r_osc(ii)^7))*(1-...
        (14*R_osc(3,ii)^2)/r_osc(ii)^2 + (21*R_osc(3,ii)^4)/r_osc(ii)^4);
    a_J4(3,ii) = (15*J4*mu_e*r_e^4*R_osc(3,ii)/(8*r_osc(ii)^7))*(5-...
        (70*R_osc(3,ii)^2)/(3*r_osc(ii)^2) + (21*R_osc(3,ii)^4)/r_osc(ii)^4);

    %J5 perterbations
    a_J5(1,ii) = (3*J5*mu_e*r_e^5*R_osc(1,ii)*R_osc(3,ii)/(8*r_osc(ii)^9))*(35-...
        210*R_osc(3,ii)^2/r_osc(ii)^2 + 231*R_osc(3,ii)^4/r_osc(ii)^4);
    a_J5(2,ii) = (3*J5*mu_e*r_e^5*R_osc(2,ii)*R_osc(3,ii)/(8*r_osc(ii)^9))*(35-...
        210*R_osc(3,ii)^2/r_osc(ii)^2 + 231*R_osc(3,ii)^4/r_osc(ii)^4);
    a_J5(3,ii) = (3*J5*mu_e*r_e^5*R_osc(3,ii)^2/(8*r_osc(ii)^9))*(35-...
        210*R_osc(3,ii)^2/r_osc(ii)^2 + 231*R_osc(3,ii)^4/r_osc(ii)^4) -...
        15*J5*mu_e*r_e^5/(8*r_osc(ii)^7);

    %J6 perterbations
    a_J6(1,ii) = (-J6*mu_e*r_e^6*R_osc(1,ii)/(16*r_osc(ii)^9))*(35-...
        945*R_osc(3,ii)/r_osc(ii)^2 + 3465*R_osc(3,ii)^4/r_osc(ii)^4 -...
        3003*R_osc(3,ii)^6/r_osc(ii)^6);
    a_J6(2,ii) = (-J6*mu_e*r_e^6*R_osc(2,ii)/(16*r_osc(ii)^9))*(35-...
        945*R_osc(3,ii)/r_osc(ii)^2 + 3465*R_osc(3,ii)^4/r_osc(ii)^4 -...
        3003*R_osc(3,ii)^6/r_osc(ii)^6);
    a_J6(3,ii) = (-J6*mu_e*r_e^6*R_osc(3,ii)/(16*r_osc(ii)^9))*(245-...
        2205*R_osc(3,ii)/r_osc(ii)^2 + 4851*R_osc(3,ii)^4/r_osc(ii)^4 -...
        3003*R_osc(3,ii)^6/r_osc(ii)^6);

        %N-Body perterbations 
    a_n(:,ii) = ((mu_s*((Fs*rsun)-R(:,ii)))/(norm(Rss)^3) + ...
        (mu_J*((FJ*rJ)-R(:,ii)))/(norm(RssJ)^3) + ...
        (mu_V*((FV*rV)-R(:,ii)))/(norm(RssV)^3));
    
        %SRP perterbations
    sdw = shadowfunction(JDf,rsun,R_osc(:,ii));
    a_srp(:,ii) = SolarRadiationPressure(R(:,ii),rsun,sdw,A_sc,m);
    
        %Total Perterbations
    ap(:,ii) = a_n(:,ii)*onoff(1) + ...
        a_drag(:,ii)*onoff(2) +...
        a_J2(:,ii)*onoff(3) + ...
        a_J3(:,ii)*onoff(4) + ...
        a_J4(:,ii)*onoff(5) + ...
        a_J5(:,ii)*onoff(6) + ...
        a_J6(:,ii)*onoff(7) + ...
        a_srp(:,ii)*onoff(8);
 
    q = dot(dr(:,ii),((2*R(:,ii))-dr(:,ii)))/(r(ii)^2) ; 
    F = (((q^2) - (3*q) + 3)/(1 + (1 - q)^(3/2))) * q ; 

    da(:,ii+1) = ((-mu_e*(dr(:,ii) - (F*R(:,ii))))/(r_osc(ii)^3)) + ap(:,ii) ;
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

    for mm = 1:length(R(1,:)) 
      [h(mm), inc(mm), RAAN(mm), arg(mm), ecc(mm), ra(mm), rp(mm)] = state2coes(R(1:3,mm), V(1:3,mm), mu_e) ;
    end
    dcoes = [h', inc', RAAN', arg', ecc', ra', rp'] ;
    
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

%Curtis' Standard Atmosphere Exponential Interpolation
function rho = atmosphere(z)
%inputs: altitude (km)
%outputs: density (kg/m^3)

%Altitudes (km)
h = [ 0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 ...
    300 350 400 450 500 600 700 800 900 1000 ];

%Corresponding Densities (kg/m^3)
r = [ 1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 ...
    1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 ...
    2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 ...
    1.184e-12 5.215e-13 1.137e-13 3.07e-14 1.136e-14 5.759e-15 3.561e-15 ];

%Scale Height (km)
H = [ 7.31 6.427 6.546 7.36 8.342 7.583 6.661 5.927 5.533 5.703 6.782 ...
    9.973 13.243 16.322 21.652 27.974 34.934 43.342 49.755 54.513 ...
    58.019 60.98 65.654 76.377 100.587 147.203 208.02 ];

%If altitude is outside of the range
if z > 1000
    z = 1000;
elseif z < 0
    z = 0;
end

%Interpolation Interval
for j = 1:27
    if z >= h(j) && z < h(j+1)
        i = j;
    end
end
if z >= 1000
    i = 27;
end

%Exponential Interpolation
rho = r(i)*exp(-(z-h(i))/H(i));

end

%Solar Radiation Pressure
function a_srp = SolarRadiationPressure(rsc,rsun,F,Asc,m)
%inputs: radius from earth to sc, radius from earth to sun, shadow
%       function, Area exposed to the sun (m^2)
%outputs: acceleration due to SRP (km/s)

%assumptions
Psr = 4.57e-6;      %N/m^2
Cr = 1.2;

rrel = rsun - rsc;

a_srp = -Psr*Cr*Asc/m * rrel/norm(rrel) * F * 10^-3;

end

%Shadow Function
function F = shadowfunction(~,rsun,rsc)
%inputs: Julian Date, position of spacecraft (ECI)
%outputs: Shadow Function for SRP Effects Calculation

    Rsun = norm(rsun);
    Rsc = norm(rsc);
    r_e = 6378;     %km
    
    theta = acosd((dot(rsc,rsun))/(Rsc*Rsun));
    theta1 = acosd(r_e/Rsc);
    theta2 = acosd(r_e/Rsun);
    thetatot = theta1 + theta2;
    
    if thetatot <= theta  
        F = 0;
    else
        F = 1;
    end
end

function [xx,yy,zz] = earth_sphere(varargin)
%EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 
%% Input Handling
[cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
% Handle remaining inputs.
% Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
for i = 1:nargs
    if ischar(args{i})
        units = args{i};
        j = j+1;
    elseif isnumeric(args{i})
        n = args{i};
        k = k+1;
    end
end
if j > 1 || k > 1
    error('Invalid input types')
end
%% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = 6378.1363*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end
     

% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
%% Plotting
if nargout == 0
    cax = newplot(cax);
    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; %#ok<NODEF>
    

    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    % Create the sphere with Earth topography and adjust colormap
    surface(x,y,z,props,'parent',cax)
    colormap(topomap1)
% Replace the calls to surface and colormap with these lines if you do 
% not want the Earth's topography displayed.
%     surf(x,y,z,'parent',cax)
%     shading flat
%     colormap gray
    

    % Refine figure
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])
    view(127.5,30)
else
    xx = x; yy = y; zz = z;
end
end
