%%
% - Rebecca Gross -
% - Aero 452 -
% - Homework 3 -
% - 11/13/19 -
%%
clear all; close all; clc;
%% Variables
r_e = 6378;
mu = 398600;
%%
disp('Problem 2')

zp0 = 300;      %perigee altitude (km)
za0 = 3092;     %apogee altitude (km)
raan0 = 45;     %right ascention of ascending node (deg)
inc0 = 28;      %inclination (deg)
omega0 = 30;    %argument of perigee (deg)
theta0 = 40;    %true anomoly (deg)

rp0 = r_e + zp0;        %radius of perigee (km)
ra0 = r_e + za0;        %radius of apogee (km)
ecc0 = (ra0-rp0)/(ra0+rp0);     %eccentricity
a0 = (ra0+rp0)/2;             %semimajor axis (km)
h0 = sqrt(a0*mu*(1-ecc0^2));    %angular momentum
T0 = 2*pi/sqrt(mu) * a0^1.5;    %period (seconds)

[r0,v0] = COES_RV(ecc0,h0,theta0,omega0,raan0,inc0);

% ------------------------ Cowell's Only J2 -------------------------------
disp(' ')
tspan = [0 48*60*60];
state = [r0',v0'];
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
onoff = [1 0];      %(1)J2     (2)J3
[tnew,statenew] = ode45(@CowellsMethod,tspan,state,options,mu,onoff);

rnew = statenew(:,1:3);
vnew = statenew(:,4:6);
h = zeros(size(tnew));
inc = zeros(size(tnew));
raan = zeros(size(tnew));
ecc = zeros(size(tnew));
omega = zeros(size(tnew));
theta = zeros(size(tnew));
rp = zeros(size(tnew));
ra = zeros(size(tnew));
for i = 1:length(tnew)
    
    [h1,inc1,raan1,ecc1,omega1,theta1,rp1,ra1] = ...
        RV_COES(rnew(i,:),vnew(i,:));
    
    h(i) = h1;
    inc(i) = inc1;
    raan(i) = raan1;
    ecc(i) = ecc1;
    omega(i) = omega1;
    theta(i) = theta1;
    rp(i) = rp1;
    ra(i) = ra1;
    
end

figure(1)
title('Cowell''s Method')
subplot(3,1,1)
plot(tnew/3600,raan)
xlabel('time (hours)')
ylabel('RAAN (degrees)')
grid on
title('RAAN vs Time')
subplot(3,1,2)
plot(tnew/3600,inc)
xlabel('time (hours)')
ylabel('Inclination (degrees)')
grid on
title('Inclination vs Time')
subplot(3,1,3)
plot(tnew/3600,omega)
xlabel('time (hours)')
ylabel('Argument of Perigee (degrees)')
grid on
title('Argument of Perigee vs Time')

figure(2)
plot(tnew/3600,rp-r_e)
hold on
plot(tnew/3600,ra-r_e)
xlabel('time (hours)')
ylabel('Altitude (km)')
grid on
title('Apogee and Perigee Altitudes')
legend('Perigee Altitude','Apogee Altitudes')

% -------------------------- Cowell's J2 and J3 ---------------------------
onoff = [1 1];      %(1)J2     (2)J3
[tnew2,statenew2] = ode45(@CowellsMethod,tspan,state,options,mu,onoff);
disp(' ')

rnew2 = statenew2(:,1:3);
vnew2 = statenew2(:,4:6);
h2 = zeros(size(tnew2));
inc2 = zeros(size(tnew2));
raan2 = zeros(size(tnew2));
ecc2 = zeros(size(tnew2));
omega2 = zeros(size(tnew2));
theta2 = zeros(size(tnew2));
rp2 = zeros(size(tnew2));
ra2 = zeros(size(tnew2));
for i = 1:length(tnew2)
    
    [h12,inc12,raan12,ecc12,omega12,theta12,rp12,ra12] = ...
        RV_COES(rnew2(i,:),vnew2(i,:));
    
    h2(i) = h12;
    inc2(i) = inc12;
    raan2(i) = raan12;
    ecc2(i) = ecc12;
    omega2(i) = omega12;
    theta2(i) = theta12;
    rp2(i) = rp12;
    ra2(i) = ra12;
    
end

figure(3)
title('Cowell''s Method')
subplot(3,1,1)
plot(tnew2/3600,raan2)
xlabel('time (hours)')
ylabel('RAAN (degrees)')
grid on
title('RAAN vs Time')
subplot(3,1,2)
plot(tnew2/3600,inc2)
xlabel('time (hours)')
ylabel('Inclination (degrees)')
grid on
title('Inclination vs Time')
subplot(3,1,3)
plot(tnew2/3600,omega2)
xlabel('time (hours)')
ylabel('Argument of Perigee (degrees)')
grid on
title('Argument of Perigee vs Time')

figure(4)
plot(tnew2/3600,rp2-r_e)
hold on
plot(tnew2/3600,ra2-r_e)
xlabel('time (hours)')
ylabel('Altitude (km)')
grid on
title('Apogee and Perigee Altitudes')
legend('Perigee Altitude','Apogee Altitudes')

disp(['There is no significant difference on the COES plots when using ',...
    'only J2 versus using both J2 and J3. This supports the idea that ',...
    'J2 is the most significant non-spherical Earth perterbation and when ',...
    'doing non-super important calculations, J3 does not need to be ',...
    'considered. NSE causes RAAN to decrease, Argument of Perigee to ',...
    'increase, and does not affect the inclination other than small oscillations.'])

%% Functions

%r and v vectors to COES
function [h,inc,RAAN,ecc,arg_p,theta,rp,ra] = RV_COES(r,v)

% Initialize constants
mu_e = 398600;

% Calculate distance and speed of S/C
R = norm(r);
V = norm(v);

% Calculate radial velocity
v_r = dot(r,v)/R;

% Calculate specific angular momentum (h)
h_hat = cross(r',v');
h = norm(h_hat);

% Calculate Inclination (inc)
inc = acosd(h_hat(3)/h);

% Calculate Right Ascension of Ascending Node (RAAN)
K_hat = [0 0 1];
N_hat = cross(K_hat,h_hat);
N = norm(N_hat);
RAAN = acosd(N_hat(1)/N);

if N_hat(2) < 0
    RAAN = 360-RAAN;
end

% Calculate eccentricity (ecc)
ecc_hat = cross(v,h_hat)/mu_e - r/R;
% ecc_hat = (1/mu_e)*(((V^2-(mu_e/R)).*r)-((R*v_r).*v));
ecc = norm(ecc_hat);

if inc ~= 0
    if ecc > eps
        arg_p = acosd(dot(N_hat,ecc_hat)/(N*ecc));
        if ecc == 0
            arg_p = 360 - w;
        end
    else
        arg_p = 0;
    end
end

% Calculate True Anomaly (theta)
if v_r >= 0
    theta = acosd(dot(ecc_hat,r)/(ecc*R));
elseif v_r < 0
    theta = 360-(dot(ecc_hat,r)/(ecc*R));
end

%Calculating radius of apogee and perigee
rp = (h^2)/(mu_e*(1+ecc));
ra = (h^2)/(mu_e*(1-ecc));

end

%COES to r and v vectors
function [r_ECI,v_ECI] = COES_RV(ecc,h,theta,arg_p,RAAN,inc)

% Initialize constants
mu_earth = 398600;   % km^3/s^2

% Find r and v vectors in PERI
r_PERI = (h^2/mu_earth)*(1/(1+ecc*cosd(theta))).* ... 
    [cosd(theta);sind(theta);0];
v_PERI = (mu_earth/h).*[-sind(theta);ecc+cosd(theta);0];

% Rotation matrices for each angle
R3a = [cosd(arg_p) sind(arg_p) 0; -sind(arg_p) cosd(arg_p) 0; 0 0 1];
R1i = [1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cos(inc)];
R3r = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];

% Find rotation matrix Q using 3-1-3 rotaion sequence
Q = R3a*R1i*R3r;

% Find r and v vectors in ECI
r_ECI = Q'*r_PERI;
v_ECI = Q'*v_PERI;

end

%Cowell's Method
function dstatedt = CowellsMethod(~,state,mue,onoff)

%%%% constants %%%%
J2 = 0.00108263;
J3 = (-2.53266e-6)*J2;
r_e = 6378;                     %radius of the earth (km)

%dr vector
dx = state(4);
dy = state(5);
dz = state(6);
r = [state(1) state(2) state(3)]';
R = norm([state(1) state(2) state(3)]);

%%%% perterbations %%%%
%J2 perterbations
temp = -((3*J2*mue*(r_e^2))/(2*(R^5)));
aI = temp*r(1)*(1 - ((5*r(3)^2)/(R^2)));
aJ = temp*r(2)*(1 - ((5*r(3)^2)/(R^2)));
aK = temp*r(3)*(3 - ((5*r(3)^2)/(R^2)));
a_J2 = [aI , aJ , aK];

%J3 perterbations
a_J3(1) = -(5*J3*mue*r_e^3*r(1)/(2*R^7))*(3*r(3) - (7*r(3)^3)/R^2);
a_J3(2) = -(5*J3*mue*r_e^3*r(2)/(2*R^7))*(3*r(3) - (7*r(3)^3)/R^2);
a_J3(3) = -(5*J3*mue*r_e^3*r(3)/(2*R^7))*(6*r(3)^2 - ((7*r(3)^4)/R^2) - (3*R^2)/5);

%total perterbations
a_pert = a_J2*onoff(1) + a_J3*onoff(2);            %km/s^2

%%%% new state vector %%%%
ddx = -mue*state(1)/R^3 + a_pert(1);
ddy = -mue*state(2)/R^3 + a_pert(2); 
ddz = -mue*state(3)/R^3 + a_pert(3);

dstatedt = [dx;dy;dz;ddx;ddy;ddz];

end
