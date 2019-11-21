%%
% - Rebecca Gross -
% - Aero 452 -
% - Homework 3 -
% - 11/2/19 -
%%
close all; clear all; clc;
%% Constants

global mu

r_earth = 6378;         %radius of the earth (km)
mu = 398600;            %km^3/s^2
% %% Problem 1
% disp('Problem 1')
% 
 t1 = 120*24*60*60;      %120 days in seconds
diam = 1;               %diameter of spacecraft (m)
m = 100;                %mass (kg)
omega_e = [0 0 72.9211e-6];     %angular velocity of the earth (rad/s)
A_sc = pi*(diam/2)^2;                   %spacecraft area (circular) (m^2)
CD = 2.2;                               %coefficient of drag

altp = 215;                 %perigee altitude (km)
rp = r_earth + altp;        %perigee radius (km)
alta = 939;                 %apogee altitude (km)
ra = r_earth + alta;        %apogee radius (km)
raan0 = deg2rad(340);       %radians
inc0 = deg2rad(65.2);       %radians
omega0 = deg2rad(58);       %radians
theta0 = deg2rad(332);      %radians
ecc0 = (ra-rp)/(ra+rp);     %eccentricity
a0 = (ra+rp)/2;             %semimajor axis (km)
h0 = sqrt(a0*mu*(1-ecc0^2));    %angular momentum
T0 = 2*pi/sqrt(mu) * a0^1.5;    %period (seconds)

 [r0,v0] = COES_RV(ecc0,h0,theta0,omega0,raan0,inc0);
 tspan = [0 t1];
 
%% ---------------------------- Cowell's ----------------------------------
disp(' ')
disp('Running Cowell''s')
state_c = [r0',v0'];
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
tic;
[tnew_c,statenew_c] = ode45(@CowellsMethod,tspan,state_c,options,mu);
toc;
t_Cowells = toc;
disp(' ')

disp('Running RV_COES')
tic;
rnew_c = statenew_c(:,1:3);
vnew_c = statenew_c(:,4:6);
h_c = zeros(size(tnew_c));
inc_c = zeros(size(tnew_c));
raan_c = zeros(size(tnew_c));
ecc_c = zeros(size(tnew_c));
omega_c = zeros(size(tnew_c));
theta_c = zeros(size(tnew_c));
rp_c = zeros(size(tnew_c));
ra_c = zeros(size(tnew_c));
for i = 1:length(tnew_c)
    
    [h,inc,raan,ecc,omega,theta,rp,ra] = ...
        RV_COES(rnew_c(i,:),vnew_c(i,:));
    
    h_c(i) = h;
    inc_c(i) = inc;
    raan_c(i) = raan;
    ecc_c(i) = ecc;
    omega_c(i) = omega;
    theta_c(i) = theta;
    rp_c(i) = rp;
    ra_c(i) = ra;
    
end
toc;
t_C_coes = toc;
disp(' ')

ind_c = find(ecc_c<0.001,1);

figure(1)
plot3(statenew_c(1:ind_c,1),statenew_c(1:ind_c,2),statenew_c(1:ind_c,3))
grid on
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Propogation of Orbit')
grid on

figure(2)
title('Cowell''s Method')
subplot(3,1,1)
plot(tnew_c(1:ind_c)/86400,raan_c(1:ind_c))
xlabel('time (days)')
ylabel('RAAN (degrees)')
grid on
title('RAAN vs Time')
subplot(3,1,2)
plot(tnew_c(1:ind_c)/86400,inc_c(1:ind_c))
xlabel('time (days)')
ylabel('Inclination (degrees)')
grid on
title('Inclination vs Time')
subplot(3,1,3)
plot(tnew_c(1:ind_c)/86400,omega_c(1:ind_c))
xlabel('time (days)')
ylabel('Argument of Perigee (degrees)')
grid on
title('Argument of Perigee vs Time')

figure(3)
plot(tnew_c(1:ind_c)/86400,rp_c(1:ind_c))
hold on
plot(tnew_c(1:ind_c)/86400,ra_c(1:ind_c))
xlabel('time (days)')
ylabel('radius (km)')
grid on
title('Radius of Apogee and Perigee')
legend('Radius of Perigee','Radius of Apogee')

t_c_tot = t_Cowells + t_C_coes;

%% ----------------------------- Encke's -----------------------------------
disp('Running Encke''s')
tic;
%universal variable crap
dt = T0/100;          %time step for UV function (sec)
time(1) = 0;           %initial time
i = 1;              %counter

%initial osc vectors
r_t = r0';
v_t = v0';
R_t = norm(r_t);
altell = 600;

[he,ince,raane,ecce,omegae,thetae,rpe,rae] = RV_COES(r_t,v_t);

while time(i) < t1 && altell > 120
    
    dr0 = [0 0 0]';
    dv0 = [0 0 0]';
    dR = norm(dr0);
    dV = norm(dv0);
    
    r_osc = r_t(i,:);
    v_osc = v_t(i,:);
    R_osc = norm(r_osc);
    
    q = 0;
    F_q = 0;
    
    R_osc = norm(r_osc);
    V_osc = norm(v_osc);
    
    %Drag perts
    vrel = (v_osc' - cross(omega_e,r_osc)');   %relative velocity (km/s)
    A_sc = (pi*(diam/2)^2)/(1000^2);        %spacecraft area (circular) (km^2)
    CD = 2.2;                               %coefficient of drag

    altell = (R_osc-6378);        %elliptical altitude (m)
    rho = (atmosphere(altell))*10^9;       %kg/km^3

    a_drag = -(1/2)*CD*(A_sc/m)*rho*(norm(vrel)^2).*vrel/norm(vrel); %ECI (km/s^2)

    a_pert = a_drag';            %km/s^2
    
    b = (mu/R_osc^3).*(dr0' - F_q.*r_t(i,:));
    
    da = a_pert - b;
    dv = da'*dt + dv0;
    dr = 0.5*da'*dt^2 + dv*dt + dr0;

    [r_osc,v_osc] = UV_Propogat(r_osc,v_osc,dt);
    
    R_osc = norm(r_osc);
    V_osc = norm(v_osc);
        
    r_t(i+1,:) = r_osc + dr';
    v_t(i+1,:) = v_osc + dv';
        
    [he(i+1),ince(i+1),raane(i+1),ecce(i+1),omegae(i+1),thetae(i+1),...
        rpe(i+1),rae(i+1)] = RV_COES(r_t(i+1,:),v_t(i+1,:));
    a0 = (rae(i+1)+rpe(i+1))/2;
    
    %setting constraints for loop
    time(i+1) = time(i)+dt;
    i = i+1;
    
end
toc;
t_Enckes = toc;
disp(' ')

figure(4)
subplot(3,1,1)
plot(time/86400,rad2deg(raane))
xlabel('time (days)')
ylabel('RAAN (degrees)')
grid on
title('Encke''s: RAAN vs Time')
subplot(3,1,2)
plot(time/86400,rad2deg(ince))
xlabel('time (days)')
ylabel('Inclination (degrees)')
grid on
title('Inclination vs Time')
subplot(3,1,3)
plot(time/86400,rad2deg(omegae))
xlabel('time (days)')
ylabel('Argument of Perigee (degrees)')
grid on
title('Argument of Perigee vs Time')

figure(5)
plot(time/86400,rpe)
hold on
plot(time/86400,rae)
xlabel('time (days)')
ylabel('radius (km)')
grid on
title('Encke''s: Radius of Apogee and Perigee')
legend('Radius of Perigee','Radius of Apogee')

%% --------------------------------- VoP -----------------------------------
disp('Running VoP')
tic
tspan = [0 120*24*60*60];
state = [ecc0,h0,theta0,omega0,raan0,inc0];
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Event',@stooopp);  
[tnew_v,statenew_v] = ode45(@GaussianVoP,tspan,state,options,mu);
toc;
t_VoP = toc;
disp(' ')

for i = 1:length(tnew_v)
    
    %Calculating radius of apogee and perigee
    rp_v(i) = (statenew_v(i,2)^2)/(mu*(1+statenew_v(i,1)));
    ra_v(i) = (statenew_v(i,2)^2)/(mu*(1-statenew_v(i,1)));
    
end

figure(6)
subplot(3,1,1)
plot(tnew_v/86400,rad2deg(statenew_v(:,5)))
xlabel('time (days)')
ylabel('RAAN (degrees)')
grid on
title('VoP: RAAN vs Time')
subplot(3,1,2)
plot(tnew_v/86400,rad2deg(statenew_v(:,6)))
xlabel('time (days)')
ylabel('Inclination (degrees)')
grid on
title('Inclination vs Time')
subplot(3,1,3)
plot(tnew_v/86400,rad2deg(statenew_v(:,4)))
xlabel('time (days)')
ylabel('Argument of Perigee (degrees)')
grid on
title('Argument of Perigee vs Time')

figure(7)
plot(tnew_v/86400,rp_v)
hold on
plot(tnew_v/86400,ra_v)
xlabel('time (days)')
ylabel('radius (km)')
grid on
title('VoP: Radius of Apogee and Perigee')
legend('Radius of Perigee','Radius of Apogee')

%% -------------------------------- 2body ----------------------------------

state_2 = [r0',v0'];
disp('Running 2-Body')
state_c = [r0',v0'];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
tic;
[tnew_2,statenew_2] = ode45(@Aero351twobodymotion,tspan,state_c,options,mu);
toc;
t_2body = toc;
disp(' ')

rnew_2 = statenew_2(:,1:3);
vnew_2 = statenew_2(:,4:6);
h_2 = zeros(size(tnew_2));
inc_2 = zeros(size(tnew_2));
raan_2 = zeros(size(tnew_2));
ecc_2 = zeros(size(tnew_2));
omega_2 = zeros(size(tnew_2));
theta_2 = zeros(size(tnew_2));
rp_2 = zeros(size(tnew_2));
ra_2 = zeros(size(tnew_2));
for i = 1:length(tnew_2)
    
    [h,inc,raan,ecc,omega,theta,rp,ra] = ...
        RV_COES(rnew_2(i,:),vnew_2(i,:));
    
    h_2(i) = h;
    inc_2(i) = inc;
    raan_2(i) = raan;
    ecc_2(i) = ecc;
    omega_2(i) = omega;
    theta_2(i) = theta;
    rp_2(i) = rp;
    ra_2(i) = ra;
    
end

figure(8)
title('2-Body Method')
subplot(3,1,1)
plot(tnew_2/86400,raan_2)
xlabel('time (days)')
ylabel('RAAN (degrees)')
ylim([5.9 6.0])
grid on
title('RAAN vs Time')
subplot(3,1,2)
plot(tnew_2/86400,inc_2)
xlabel('time (days)')
ylabel('Inclination (degrees)')
ylim([1.13 1.14])
grid on
title('Inclination vs Time')
subplot(3,1,3)
plot(tnew_2/86400,omega_2)
xlabel('time (days)')
ylabel('Argument of Perigee (degrees)')
ylim([1 1.1])
grid on
title('Argument of Perigee vs Time')

figure(9)
plot(tnew_2/86400,rp_2)
hold on
plot(tnew_2/86400,ra_2)
xlabel('time (days)')
ylabel('radius (km)')
grid on
title('2-Body: Radius of Apogee and Perigee')
legend('Radius of Perigee','Radius of Apogee')

disp(['The total time Cowell''s Method took to run was ',num2str(t_c_tot),...
    ' seconds. Encke''s took ',num2str(t_Enckes),' seconds. VoP took ',...
    num2str(t_VoP),' seconds. 2-Body took ',num2str(t_2body),' seconds.',...
    ' When comparing these times, we see that Cowell''s took the longest, ',...
    'followed by VoP, then Encke''s, and finally 2-body. All methods ',...
    'result in similar graphs and trends.'])

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
function dstatedt = CowellsMethod(~,state,mue)

%%%% constants %%%%
diam = 1;                       %diameter of spacecraft (m)
m = 100;                        %mass (kg)
omega_e = [0 0 72.9211e-6]';    %angular velocity of the earth (rad/s)

%dr vector
dx = state(4);
dy = state(5);
dz = state(6);
dr = [state(4) state(5) state(6)]';
r = [state(1) state(2) state(3)]';
R = norm([state(1) state(2) state(3)]);

%%%% perterbations %%%%
%Drag perts
vrel = (dr - cross(omega_e,r))*10^3;   %relative velocity (m/s)
A_sc = pi*(diam/2)^2;                   %spacecraft area (circular) (m^2)
CD = 2.2;                               %coefficient of drag

altell = (R-6378);        %elliptical altitude (m)
rho = atmosphere(altell);       %kg/m^3

a_drag = -(1/2)*CD*(A_sc/m)*rho*(norm(vrel)^2)*vrel/norm(vrel); %ECI (m/s^2)
a_drag = a_drag*10^-3;      %km/s^2

a_pert = a_drag;            %km/s^2

%%%% new state vector %%%%
ddx = -mue*state(1)/R^3 + a_pert(1);
ddy = -mue*state(2)/R^3 + a_pert(2); 
ddz = -mue*state(3)/R^3 + a_pert(3);

dstatedt = [dx;dy;dz;ddx;ddy;ddz];

end

function [R,V] = UV_Propogat(R0, V0, t)
    %CREDIT -- CURTIS
 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %{
      This function computes the state vector (r,v) from the 
      initial state vector (r0,v0) and the change in true anomaly.
 

      mu - gravitational parameter (km^3/s^2)
      r0 - initial position vector (km)
      v0 - initial velocity vector (km/s)
      dt - time since x = 0 (s)
      r  - final position vector (km)
      v  - final velocity vector (km/s)
 

      User M-functions required: f_and_g_ta, fDot_and_gDot_ta
    %}
    % --------------------------------------------------------------------
 

    global mu
 

    %...Magnitudes of R0 and V0: 
        r0 = norm(R0);
        v0 = norm(V0);
    %...Initial radial velocity: 
        vr0 = dot(R0, V0)/r0;
    %...Reciprocal of the semimajor axis (from the energy equation):   
        alpha = 2/r0 - v0^2/mu;
    %...Compute the universal anomaly: 
        x = kepler_U(t, r0, vr0, alpha);
 

 

 

    %...Compute the f and g functions and their derivatives:
        [f, g]       =       f_and_g(x, t, r0, alpha);
    %...Compute the final position vector: 
        R = f*R0 + g*V0;
    %...Compute the magnitude of R: 
        r = norm(R);
    %...Compute the derivatives of f and g:
        [fdot, gdot] = fDot_and_gDot(x, r, r0, alpha);
 

    %...Compute the final position and velocity vectors:
        V = fdot*R0 + gdot*V0;
 

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function [f, g] = f_and_g(x, t, ro, a)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %{
      This function calculates the Lagrange f and g coefficients.
 

      mu - the gravitational parameter (km^3/s^2)
      a  - reciprocal of the semimajor axis (1/km)
      ro - the radial position at time to (km)
      t  - the time elapsed since ro (s)
      x  - the universal anomaly after time t (km^0.5)
      f  - the Lagrange f coefficient (dimensionless)
      g  - the Lagrange g coefficient (s)
 

      User M-functions required:  stumpC, stumpS
    %}
    % ----------------------------------------------
 

    global mu
 

    z = a*x^2;
 

    %...Equation 3.69a:
    f = 1 - x^2/ro*stumpC(z);
 

    %...Equation 3.69b:
    g = t - 1/sqrt(mu)*x^3*stumpS(z);
 

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function [fdot, gdot] = fDot_and_gDot(x, r, ro, a)
    %CREDIT -- CURTIS
 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %{
      This function calculates the time derivatives of the
      Lagrange f and g coefficients.
 

      mu    - the gravitational parameter (km^3/s^2)
      a     - reciprocal of the semimajor axis (1/km)
      ro    - the radial position at time to (km)
      t     - the time elapsed since initial state vector (s)
      r     - the radial position after time t (km)
      x     - the universal anomaly after time t (km^0.5)
      fdot  - time derivative of the Lagrange f coefficient (1/s)
      gdot  - time derivative of the Lagrange g coefficient (dimensionless)
 

      User M-functions required:  stumpC, stumpS
    %}
    % --------------------------------------------------
 

    global mu
 

    z = a*x^2;
 

    %...Equation 3.69c:
    fdot = sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;
 

    %...Equation 3.69d:
    gdot = 1 - x^2/r*stumpC(z);
 

    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function x = kepler_U(dt, ro, vro, a)
    %CREDIT -- CURTIS
 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %{
 

    %Inputs:
    %  mu  - gravitational parameter (km^3/s^2)
       x   - the universal anomaly (km^0.5)
       dt  - time since x = 0 (s)
       ro  - radial position when x = 0 (km)
       vro - radial velocity when x = 0 (km/s)
       a   - semimajor axis (km)
 

 

      This function uses Newton's method to solve the universal
      Kepler equation for the universal anomaly.
 

      mu   - gravitational parameter (km^3/s^2)
      x    - the universal anomaly (km^0.5)
      dt   - time since x = 0 (s)
      ro   - radial position (km) when x = 0
      vro  - radial velocity (km/s) when x = 0
      a    - reciprocal of the semimajor axis (1/km)
      z    - auxiliary variable (z = a*x^2)
      C    - value of Stumpff function C(z)
      S    - value of Stumpff function S(z)
      n    - number of iterations for convergence
      nMax - maximum allowable number of iterations
 

      User M-functions required: stumpC, stumpS
    %}
    % ----------------------------------------------
    global mu
 

    %...Set an error tolerance and a limit on the number of iterations:
    error = 1.e-8;
    nMax  = 1000;
 

    %...Starting value for x:
    x = sqrt(mu)*abs(a)*dt;
 

    %...Iterate on Equation 3.65 until until convergence occurs within
    %...the error tolerance:
    n     = 0;
    ratio = 1;
    while abs(ratio) > error && n <= nMax
        n     = n + 1;
        C     = stumpC(a*x^2);
        S     = stumpS(a*x^2);
        F     = ro*vro/sqrt(mu)*x^2*C + (1 - a*ro)*x^3*S + ro*x - sqrt(mu)*dt;
        dFdx  = ro*vro/sqrt(mu)*x*(1 - a*x^2*S) + (1 - a*ro)*x^2*C + ro;
        ratio = F/dFdx;
        x     = x - ratio;
    end
 

    %...Deliver a value for x, but report that nMax was reached:
    if n > nMax
        fprintf('\n **No. iterations of Kepler''s equation = %g', n)
        fprintf('\n   F/dFdx                              = %g\n', F/dFdx)
    end
 

    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function s = stumpS(z)
    %CREDIT -- CURTIS
 

    % $$$$$$$$$$$$$$$$$ %{
    % This function evaluates the Stumpff function S(z) according to Equation 3.52.
    % z - input argument s - value of S(z)
    % User M-functions required: none %}
       % ??????????????????????????????????????????????
    if z > 0
    s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
    s = 1/6;
    end
 

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function c = stumpC(z)
    % $$$$$$$$$$$$$$$$$ %{
    % This function evaluates the Stumpff function C(z) according to Equation 3.53.
    % z - input argument c - value of C(z)
    % User M-functions required: none %}
       
    if z > 0
    c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);
    else
    c = 1/2;
    end
 

 end

%VoP
function [dstate] = GaussianVoP(t,state,mu)
%inputs: state vector
%outputs: derivative of state vector

%constants
diam = 1;                       %diameter of spacecraft (m)
m = 100;                        %mass (kg)
omega_e = [0 0 72.9211e-6];     %angular velocity of the earth (rad/s)

ecc = state(1);
h_norm = state(2);
theta = state(3);
omega = state(4);
raan = state(5);
inc = state(6);

%finding the position and velocity vectors from the COES
[r,v] = COES_RV(ecc,h_norm,rad2deg(theta),rad2deg(omega),rad2deg(raan),rad2deg(inc));

r_norm = norm(r);
h = cross(r,v);
h_norm = norm(h);

% --------------------------- Perterbations -------------------------------
%Drag Perts
vrel = (v - cross(omega_e,r)')*10^3;       %relative velocity (m/s)
A_sc = pi*(diam/2)^2;               %spacecraft area (citcular) (m^2)
CD = 2.2;                           %coefficient of drag

altell = (r_norm-6378);         %km
rho = atmosphere(altell);       %kg/m^3

a_drag = -(1/2)*CD*(A_sc/m)*rho*(norm(vrel)^2).*vrel/norm(vrel); %ECI (m/s^2)
a_drag = a_drag*10^-3;  %km/s^2

a_pert = a_drag';

Rhat = r/r_norm;
Nhat = h/h_norm;
That = cross(Nhat,Rhat)/norm(cross(Nhat,Rhat));

R = dot(a_pert,Rhat);
N = dot(a_pert,Nhat);
T = dot(a_pert,That);

% --------------------------- COEs propogation ----------------------------
    u = theta+omega;
    %Angular Momentum (h)
       dh = r_norm*T;
    %Eccentricity
       decc = (h_norm/mu)*R*sin(theta) + (T/(h_norm*mu))*((h_norm^2 + mu*r_norm)*...
           cos(theta) + mu*ecc*r_norm);
    %True Anomaly
       dtheta = (h_norm/r_norm^2) + (1/(ecc*h_norm))*((h_norm^2/mu)*cos(theta)...
           *R - (r_norm + (h_norm^2/mu))*sin(theta)*T);
    %Inclination
       dinc = (r_norm/h_norm)*cos(u)*N;
    %RAAN
       draan = (r_norm*sin(u)*N)/(h_norm*sin(inc));
    %Argument of perigree
       domega = (-1/(ecc*h_norm))*((h_norm^2/mu)*cos(theta)*R - ...
           (r_norm + h_norm^2/mu)*sin(theta)*T) ...
                   - (r_norm*N*sin(u))/(h_norm*tan(inc));

dstate = [decc,dh,dtheta,domega,draan,dinc]';

end

function [check,Terminator,direction] = stooopp(t,state,mu)
%same inputs as ODE45 VoP function
%will stop ODE if deorbits (altitude < 100 km)

%unpacking state vector
ecc = state(1);
h_norm = state(2);
theta = state(3);
omega = state(4);
raan = state(5);
inc = state(6);

%finding r and v vectors
[r,~] = COES_RV(ecc,h_norm,rad2deg(theta),rad2deg(omega),rad2deg(raan),rad2deg(inc));

alt = norm(r) - 6378;   %altitude

if alt >= 100
    check = 1;
elseif alt < 100
        check = 0;
end

Terminator = 1;

direction = 0;

end

%2-Body
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

%Curtis' Standard Atmosphere Exponential Interpolation
function rho = atmosphere(z)
%inputs: altitude (km)
%outputs: density (kg/m^3)

%Altitudes (km)
h = [ 0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 ...
    300 350 400 450 500 600 700 800 1000 ];

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
if z == 1000
    i = 27;
end

%Exponential Interpolation
rho = r(i)*exp(-(z-h(i))/H(i));

end
