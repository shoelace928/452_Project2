%% Notes to self

%Ask Dr. A:
% Method for delta-v calcs - does this seem like a good aapproach?

% Getting high numbers for delta-v burns? (~28 km/s each correction)
% Burn values are varying and that doesn't make sense?
% Good time step for corrections? (day vs period)

% Part c - is it fine to just do a delta-v calc for an inc/raan change and 
% talk about how you don't need to do that burn for sun synchronous orbits 
% or is that not enough?

%% Part 1 - Propagation

clear all; close all; clc; 
mu_e = 398600 ; 
W = [0 0 72.9211e-6]';  %angular velocity of the earth (rad/s)
r_e = 6378;             %radius of the Earth (km)
tspan = 7*24*3600 ; %sec
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

%% Iridium - Drag, N-Body, Oblateness, SRP
RI = [-4351.6; -5661.3; 15.3] ;
VI = [.3814; -.2730; 7.4567] ;
stateI = [RI VI] ;
A_i = 19.25;   %m^2
m_i = 860;  %kg
TI = 6.004910595e3 ;   %period, seconds

%Only Sun N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 0 0 1 1 1 1 1 1 1 0];
dt = 100 ; %time step in seconds 
[RI2, VI2, tt, dcoesI, ~] = enckesnb(RI, RI, VI, VI, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_i, m_i, r_e, TI)  ;

%Sun + Jupiter N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 1 0 1 1 1 1 1 1 1 0];
[RI2J, VI2J, tt, dcoesIJ, ~] = enckesnb(RI, RI, VI, VI, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_i, m_i, r_e, TI)  ;

%Sun + Venus N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 0 1 1 1 1 1 1 1 1 0];
[RI2V, VI2V, tt, dcoesIV, ~] = enckesnb(RI, RI, VI, VI, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_i, m_i, r_e, TI)  ;

%Assigning COES
incI = rad2deg(dcoesI(:,2));
raanI = rad2deg(dcoesI(:,3));
argI = rad2deg(dcoesI(:,4));

incIJ = rad2deg(dcoesIJ(:,2));
raanIJ = rad2deg(dcoesIJ(:,3));
argIJ = rad2deg(dcoesIJ(:,4));

incIV = rad2deg(dcoesIV(:,2));
raanIV = rad2deg(dcoesIV(:,3));
argIV = rad2deg(dcoesIV(:,4));

dincI(1) = 0; dincIJ(1) = 0; dincIV(1) = 0;
draanI(1) = 0; draanIJ(1) = 0; draanIV(1) = 0;
dargI(1) = 0; dargIJ(1) = 0; dargIV(1) = 0;
for jj = 2:length(tt)
    
    dincI(jj) = incI(jj) - incI(1);
    draanI(jj) = raanI(jj) - raanI(1);
    dargI(jj) = argI(jj) - argI(1);
    
    dincIJ(jj) = incIJ(jj) - incIJ(1);
    draanIJ(jj) = raanIJ(jj) - raanIJ(1);
    dargIJ(jj) = argIJ(jj) - argIJ(1);
    
    dincIV(jj) = incIV(jj) - incIV(1);
    draanIV(jj) = raanIV(jj) - raanIV(1);
    dargIV(jj) = argIV(jj) - argIV(1);
    
end

% aI = (dcoesI(1,6) + dcoesI(1,7))/2;
% TI = 2*pi/sqrt(mu_e) * aI^1.5;    %period (seconds)

figure(1)
earth_sphere
hold on
plot3(RI2(1,:),RI2(2,:),RI2(3,:))
grid on
title('Iridium Propagation')

figure(2)
earth_sphere
hold on
plot3(RI2J(1,:),RI2J(2,:),RI2J(3,:))
grid on
title('Iridium Propagation (Jupiter incl)')

figure(3)
earth_sphere
hold on
plot3(RI2V(1,:),RI2V(2,:),RI2V(3,:))
grid on
title('Iridium Propagation (Venus incl)')

figure(4)
title('Iridium COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincI)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanI)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargI)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

figure(5)
title('Iridium COES over time (JUPITER)')
subplot(3,1,1)
plot(tt/(24*3600),dincIJ)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanIJ)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargIJ)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

figure(6)
title('Iridium COES over time (VENUS)')
subplot(3,1,1)
plot(tt/(24*3600),dincIV)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanIV)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargIV)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

%delta-v calc:
RI2_end = [RI2(1,end);RI2(2,end);RI2(3,end)] ;
VI2_end = [VI2(1,end);VI2(2,end);VI2(3,end)] ;
[VI3, VI3_prime] = lambert(RI,RI2_end,TI) ;
burnI1 = VI3 - VI2_end ;
burnI2 = VI - VI3_prime ;
dVI = norm(burnI1) + norm(burnI2) ; 

%% COSMOS - N-body, Oblateness, SRP 
RC = [-10401; 39485; 10622] ;
VC = [-2.9744; -.7741; -.0349] ;
stateC = [RC VC] ;
A_c = 40.48;  %m^2
m_c = 1700;  %kg

%Sun N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 1 1 0 1 1 1 1 1 1 0];
dt = 100 ; %time step in seconds 
[RC2, VC2, tt, dcoesC, ~] = enckesnb(RC, RC, VC, VC, [0;0;0], [0;0;0], dt, ...
    tspan, mu_e, onoff, W, A_c, m_c, r_e, TI)  ;

%Sun + Jupiter N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 1 0 1 1 1 1 1 1 1 0];
[RC2J, VC2J, tt, dcoesCJ, ~] = enckesnb(RC, RC, VC, VC, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_c, m_c, r_e, TI)  ;

%Sun + Venus N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 0 1 1 1 1 1 1 1 1 0];
[RC2V, VC2V, tt, dcoesCV, ~] = enckesnb(RC, RC, VC, VC, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_c, m_c, r_e, TI)  ;

%Assigning COES
incC = rad2deg(dcoesC(:,2));
raanC = rad2deg(dcoesC(:,3));
argC = rad2deg(dcoesC(:,4));

incCJ = rad2deg(dcoesCJ(:,2));
raanCJ = rad2deg(dcoesCJ(:,3));
argCJ = rad2deg(dcoesCJ(:,4));

incCV = rad2deg(dcoesCV(:,2));
raanCV = rad2deg(dcoesCV(:,3));
argCV = rad2deg(dcoesCV(:,4));

dincC(1) = 0; dincCJ(1) = 0; dincCV(1) = 0;
draanC(1) = 0; draanCJ(1) = 0; draanCV(1) = 0;
dargC(1) = 0; dargCJ(1) = 0; dargCV(1) = 0;
for jj = 2:length(tt)
    
    dincC(jj) = incC(jj) - incC(1);
    draanC(jj) = raanC(jj) - raanC(1);
    dargC(jj) = argC(jj) - argC(1);
    
    dincCJ(jj) = incCJ(jj) - incCJ(1);
    draanCJ(jj) = raanCJ(jj) - raanCJ(1);
    dargCJ(jj) = argCJ(jj) - argCJ(1);
    
    dincCV(jj) = incCV(jj) - incCV(1);
    draanCV(jj) = raanCV(jj) - raanCV(1);
    dargCV(jj) = argCV(jj) - argCV(1);
    
end

aC = (dcoesC(1,6) + dcoesC(1,7))/2;
TC = 2*pi/sqrt(mu_e) * aC^1.5;    %period (seconds)

figure(7)
earth_sphere
hold on
plot3(RC2(1,:),RC2(2,:),RC2(3,:))
grid on
title('Cosmos Propagation')

figure(8)
earth_sphere
hold on
plot3(RC2J(1,:),RC2J(2,:),RC2J(3,:))
grid on
title('Cosmos Propagation')

figure(9)
earth_sphere
hold on
plot3(RC2V(1,:),RC2V(2,:),RC2V(3,:))
grid on
title('Cosmos Propagation')

figure(10)
title('Cosmos COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincC)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanC)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargC)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

figure(11)
title('Cosmos COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincCJ)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanCJ)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargCJ)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

figure(12)
title('Cosmos COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincCV)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanCV)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargCV)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

%delta-v calc:
RC2_end = [RC2(1,end);RC2(2,end);RC2(3,end)] ;
VC2_end = [VC2(1,end);VC2(2,end);VC2(3,end)] ;
[VC3, VC3_prime] = lambert(RC,RC2_end,TC) ;
burnC1 = VC3 - VC2_end ;
burnC2 = VC - VC3_prime ;
dVC = norm(burnC1) + norm(burnC2) ; 

%% RADARSAT-2 - N-Body, Oblateness, SRP, Drag
RR = [-2213.2 88.2 6819.2]';
VR = [-6.0621 3.8439 -2.0172]';

stateR = [RR VR] ;
A_r = 15.18;  %m^2
m_r = 2200;  %kg

%Sun N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 1 1 1 1 1 1 1 1 1 1 0];
dt = 100 ; %time step in seconds 
[RR2, VR2, tt, dcoesR, ~] = enckesnb(RR, RR, VR, VR, [0;0;0], [0;0;0], dt, ...
    tspan, mu_e, onoff, W, A_r, m_r, r_e, TI)  ;

%Sun + Jupiter N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 1 0 1 1 1 1 1 1 1 0];
[RR2J, VR2J, tt, dcoesRJ, ~] = enckesnb(RR, RR, VR, VR, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_r, m_r, r_e, TI)  ;

%Sun + Venus N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 0 1 1 1 1 1 1 1 1 0];
[RR2V, VR2V, tt, dcoesRV, ~] = enckesnb(RR, RR, VR, VR, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_r, m_r, r_e, TI)  ;

%Assigning COES
incR = rad2deg(dcoesR(:,2));
raanR = rad2deg(dcoesR(:,3));
argR = rad2deg(dcoesR(:,4));

incRJ = rad2deg(dcoesRJ(:,2));
raanRJ = rad2deg(dcoesRJ(:,3));
argRJ = rad2deg(dcoesRJ(:,4));

incRV = rad2deg(dcoesRV(:,2));
raanRV = rad2deg(dcoesRV(:,3));
argRV = rad2deg(dcoesRV(:,4));

dincR(1) = 0; dincRJ(1) = 0; dincRV(1) = 0;
draanR(1) = 0; draanRJ(1) = 0; draanRV(1) = 0;
dargR(1) = 0; dargRJ(1) = 0; dargRV(1) = 0;
for jj = 2:length(tt)
    
    dincR(jj) = incR(jj) - incR(1);
    draanR(jj) = raanR(jj) - raanR(1);
    dargR(jj) = argR(jj) - argR(1);
    
    dincRJ(jj) = incRJ(jj) - incRJ(1);
    draanRJ(jj) = raanRJ(jj) - raanRJ(1);
    dargRJ(jj) = argR(jj) - argRJ(1);
    
    dincRV(jj) = incRV(jj) - incRV(1);
    draanRV(jj) = raanRV(jj) - raanRV(1);
    dargRV(jj) = argRV(jj) - argRV(1);
    
end

aR = (dcoesR(1,6) + dcoesR(1,7))/2;
TR = 2*pi/sqrt(mu_e) * aR^1.5;    %period (seconds)

figure(13)
earth_sphere
hold on
plot3(RR2(1,:),RR2(2,:),RR2(3,:))
grid on
title('RadarSat-2 Propagation')

figure(14)
earth_sphere
hold on
plot3(RR2J(1,:),RR2J(2,:),RR2J(3,:))
grid on
title('RadarSat-2 Propagation')

figure(15)
earth_sphere
hold on
plot3(RR2V(1,:),RR2V(2,:),RR2V(3,:))
grid on
title('RadarSat-2 Propagation')

figure(16)
title('RadarSat-2 COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincR)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanR)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargR)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

figure(17)
title('RadarSat-2 COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincRJ)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanRJ)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargRJ)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

figure(18)
title('RadarSat-2 COES over time (SUN)')
subplot(3,1,1)
plot(tt/(24*3600),dincRV)
xlabel('time (days)')
ylabel('Inclination (deg)')
grid on
subplot(3,1,2)
plot(tt/(24*3600),draanRV)
xlabel('time (days)')
ylabel('RAAN (deg)')
subplot(3,1,3)
plot(tt/(24*3600),dargRV)
xlabel('time (days)')
ylabel('Argument of Perigee (deg)')

%delta-v calc:
RR2_end = [RR2(1,end);RR2(2,end);RR2(3,end)] ;
VR2_end = [VR2(1,end);VR2(2,end);VR2(3,end)] ;
[VR3, VR3_prime] = lambert(RR,RR2_end,TR) ;
burnR1 = VR3 - VR2_end ;
burnR2 = VR - VR3_prime ;
dVR = norm(burnR1) + norm(burnR2) ; 

%% Part 2 - Daily Corrections for Iridium Spacecraft

%Only Sun N-Body
%(1-3)N-Body(Sun-Jupiter-Venus) (4)Drag (5)J2 (6)J3 (7)J4 (8)J5 (9)J6
%(10)SRP (11)Daily Corrections
onoff = [1 0 0 1 1 1 1 1 1 1 1];
dt = 100 ; %time step in seconds 
[RI22, VI22, tt, dcoesI2, burn2] = enckesnb(RI, RI, VI, VI, [0;0;0], [0;0;0], ...
    dt, tspan, mu_e, onoff, W, A_i, m_i, r_e, TI)  ;

ind = find(burn2);

burn_2 = burn2(ind);
totburn = sum(burn_2);

%% Functions: 
    %Enckes Method With Perturbations -- N Body -- Drag -- J2-6 -- SRP
function [R, V, tt, dcoes, burn] = enckesnb(R, R_osc, V, V_osc, dro, dvo, dt, tf, ...
    mu_e, onoff, omega_e, A_sc, m, r_e,TI)  
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
    J2 = -0.00108263;
    J3 = 2.53266e-6;
    J4 = 1.619e-6;
    J5 = 2.273e-7;
    J6 = -5.407e-7; 
    
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
    a_nS(:,ii) = (mu_s*((Fs*rsun)-R(:,ii)))/(norm(Rss)^3);
    a_nJ(:,ii) = (mu_J*((FJ*rJ)-R(:,ii)))/(norm(RssJ)^3);
    a_nV(:,ii) = (mu_V*((FV*rV)-R(:,ii)))/(norm(RssV)^3);
    
        %SRP perterbations
    sdw = shadowfunction(JDf,rsun,R_osc(:,ii));
    a_srp(:,ii) = SolarRadiationPressure(R(:,ii),rsun,sdw,A_sc,m);
    
        %Total Perterbations
    ap(:,ii) = a_nS(:,ii)*onoff(1) + ...
        a_nJ(:,ii)*onoff(2) + ...
        a_nV(:,ii)*onoff(3) + ...
        a_drag(:,ii)*onoff(4) + ...
        a_J2(:,ii)*onoff(5) + ...
        a_J3(:,ii)*onoff(6) + ...
        a_J4(:,ii)*onoff(7) + ...
        a_J5(:,ii)*onoff(8) + ...
        a_J6(:,ii)*onoff(9) + ...
        a_srp(:,ii)*onoff(10);
 
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

burn(ii) = 0;
%Burn corrections for Iridium (Part 2)
if onoff(11) == 1
    
    %seconds in each period for a week
%     time = linspace(1,100)*(90*60);
    time = [1 2 3 4 5 6 7]*24*3600;     %once per day for a week

    for jj = 1:length(time)
        if tt(ii) == time(jj)
            
            %propogating to find ideal location using ODE
            options = odeset('RelTol',1e-8,'AbsTol',1e-8);
            tspan = [0 tt+3*TI/4];
            [statenewdvf] = ode45(@Aero351twobodymotion,tspan,[R(:,1)' V(:,1)'],options,mu_e);
            
            %desired R and V vectors after burn
            Rdvf = statenewdvf.y(1:3,end);
            Vdvf = statenewdvf.y(4:6,end);
                    
            %finding delta-v's for lamberts
            [v1,v2] = retrolambert(Rdvf,R_osc(:,ii),3*TI/4);
            dv1I(:,ii) = v1 - V_osc(:,ii);
            dv2I(:,ii) = Vdvf - v2;
            burn(ii) = norm(dv1I(:,ii)) + norm(dv2I(:,ii));
            
        end
        
        %reassigning R vector for post-Lamberts burn
        if tt(ii) == time(jj)+3*TI/4
            R(:,ii) = Rdvf;
            V(:,ii) = Vdvf;
        end
    end
else
    burn(ii) = 0;
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

%Lacey Davis
%Lambert's problem for orbits 
function [v1,v2] = lambert(r1,r2, delta_t_given)
%given r1 and r2, solve for v1 and v2
mu = 398600 ;
r1_mag = norm(r1)   ;
r2_mag = norm(r2)   ;
crossy = cross(r1,r2) ;
z_comp = crossy(3)   ;
    if z_comp >= 0 
        delta_theta = acos(dot(r1,r2)/(r1_mag*r2_mag)) ; 
    else
        delta_theta = (2*pi) - acos(dot(r1,r2)/(r1_mag*r2_mag)) ;
    end 
A = sin(delta_theta)*sqrt((r1_mag*r2_mag)/(1-cos(delta_theta))) ; 
    if A == 0
        disp('the constant is zero')
    end
    
%z guess bounds and initial conditions
z = 0 ;
zupper = 4*pi^2 ;
zlower = -4*pi^2 ;
ii = 1 ;
delta_t_loop = 10000000 ;
TOL = 10e-8 ;
    %Z VALUE THROUGH BISECTION METHOD
while abs(delta_t_loop - delta_t_given) > TOL 
    %STUMPFF FUNCTIONS
    if z(ii)>0 
        S = ((sqrt(z(ii))) - sin(sqrt(z(ii))))/((sqrt(z(ii)))^3) ;
        C = (1 - cos(sqrt(z(ii))))/z(ii) ;
    elseif z<0
        S = (sinh(sqrt(-z(ii))) - sqrt(-z(ii)))/((sqrt(-z(ii))^3)) ; 
        C = (cosh(sqrt(-z(ii))) - 1)/-z(ii) ;
    else
        S = 1/6 ;
        C = 1/2 ;
    end 
    
    %y, chi, delta_t_loop
    y = r1_mag + r2_mag + ((A*((z*S)-1))/sqrt(C)) ;
    chi = sqrt(y/C) ;  
    delta_t_loop = (((y/C)^(3/2)*S) + (A*sqrt(y)))/sqrt(mu) ;
    
    if delta_t_loop < delta_t_given 
        zlower = z ;
        z = (zupper+zlower)/2 ;
    else
        zupper = z ;
        z = (zupper+zlower)/2 ; 
    end   
end 

%lagrange multipliers
f = 1 - (y/r1_mag) ;
g = A*(sqrt(y/mu)) ;
g_dot = 1 - (y/r2_mag) ;
f_dot = ((f*g_dot) - 1) / g ;

%v1 and v2 
v1 = (1/g)*(r2-(f*r1)) ;
v2 = (f_dot*r1) + (g_dot*v1) ;

end 

function [v1,v2] = retrolambert(r1,r2, delta_t_given)
%given r1 and r2, solve for v1 and v2
mu = 398600 ;
r1_mag = norm(r1)   ;
r2_mag = norm(r2)   ;
crossy = cross(r1,r2) ;
z_comp = crossy(3)   ;
    if z_comp >= 0 
         delta_theta = (2*pi) - acos(dot(r1,r2)/(r1_mag*r2_mag)) ;
    else
        delta_theta = acos(dot(r1,r2)/(r1_mag*r2_mag)) ;
       
    end 
A = sin(delta_theta)*sqrt((r1_mag*r2_mag)/(1-cos(delta_theta))) ; 
    if A == 0
        disp('the constant is zero')
    end
    
%z guess bounds and initial conditions
z = 0 ;
zupper = 4*pi^2 ;
zlower = -4*pi^2 ;
ii = 1 ;
delta_t_loop = 10000000 ;
TOL = 10e-8 ;
    %Z VALUE THROUGH BISECTION METHOD
while abs(delta_t_loop - delta_t_given) > TOL 
    %STUMPFF FUNCTIONS
    if z(ii)>0 
        S = ((sqrt(z(ii))) - sin(sqrt(z(ii))))/((sqrt(z(ii)))^3) ;
        C = (1 - cos(sqrt(z(ii))))/z(ii) ;
    elseif z<0
        S = (sinh(sqrt(-z(ii))) - sqrt(-z(ii)))/((sqrt(-z(ii))^3)) ; 
        C = (cosh(sqrt(-z(ii))) - 1)/-z(ii) ;
    else
        S = 1/6 ;
        C = 1/2 ;
    end 
    
    %y, chi, delta_t_loop
    y = r1_mag + r2_mag + ((A*((z*S)-1))/sqrt(C)) ;
    chi = sqrt(y/C) ;  
    delta_t_loop = (((y/C)^(3/2)*S) + (A*sqrt(y)))/sqrt(mu) ;
    
    if delta_t_loop < delta_t_given 
        zlower = z ;
        z = (zupper+zlower)/2 ;
    else
        zupper = z ;
        z = (zupper+zlower)/2 ; 
    end   
end 

%lagrange multipliers
f = 1 - (y/r1_mag) ;
g = A*(sqrt(y/mu)) ;
g_dot = 1 - (y/r2_mag) ;
f_dot = ((f*g_dot) - 1) / g ;

%v1 and v2 
v1 = (1/g)*(r2-(f*r1)) ;
v2 = (f_dot*r1) + (g_dot*v1) ;

end 

%Linear Motion ODE
function dstatedt = Aero351twobodymotion(t,state,mue)

dx = state(4);
dy = state(5);
dz = state(6);

r = norm([state(1) state(2) state(3)]);

ddx = -mue*state(1)/r^3;
ddy = -mue*state(2)/r^3;
ddz = -mue*state(3)/r^3;

dstatedt = [dx;dy;dz;ddx;ddy;ddz];

end
