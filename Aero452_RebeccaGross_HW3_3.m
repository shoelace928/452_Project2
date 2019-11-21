%%
% - Rebecca Gross -
% - Aero 452 -
% - Homework 3 -
% - 11/13/19 -
%%
clear all; close all; clc;

%% Constants

global mu

r_earth = 6378;         %radius of the earth (km)
mu = 398600;            %km^3/s^2
%% Problem 3
disp('Problem 1')

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

disp('Running Encke''s')
tic;
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
    
    %time and date of mission
    utc = [2019 1 0 0 0 0];  
    date = datevec(seconds(time(i)));
    utc = utc + date;
    
    if utc(3) ==0
        utc(3) =1;
    end
    if utc(3) >30
        utc(3) = 1;
        utc(2) = utc(2) + 1;
    end
    
    %using eci position to find latitude, longitude, and altitude
    [lla] = eci2lla(r_osc*10^3,utc);

    lat = lla(1);               %latitude
    long = lla(2);              %longitude
    alt = lla(3);               %altitude (km)

    [~,rho_all] = atmosnrlmsise00(alt,lat,long,utc(1),utc(3),utc(6),'none');
    rho = rho_all(6) * 10^9;    %kg/km^3
    
    a_drag = -(1/2)*CD*(A_sc/m)*rho*(norm(vrel)^2).*vrel/norm(vrel); %ECI (km/s^2)

    a_pert = a_drag';            %km/s^2
    
    b = (mu/R_osc^3).*(dr0' - F_q.*r_t(i,:));
    
    da = a_pert - b;
    dv = da'*dt + dv0;
    dr = 0.5*da'*dt^2 + dv*dt + dr0;

    %Universal Variable to propogate
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

figure(1)
title('Encke''s Method')
subplot(3,1,1)
plot(time/86400,rad2deg(raane))
xlabel('time (days)')
ylabel('RAAN (degrees)')
grid on
title('RAAN vs Time')
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

figure(2)
plot(time/86400,rpe)
hold on
plot(time/86400,rae)
xlabel('time (days)')
ylabel('radius (km)')
grid on
title('Radius of Apogee and Perigee')
legend('Radius of Perigee','Radius of Apogee')

disp(['When using the higher fidelity density model, the plots are very ',...
    'similar to those produced from the exponential model. The higher ',...
    'fidelity model makes the spacecraft deorbit sooner, at about 104 ',...
    'days vs about 115 days with the exponential model. The inclination ',...
    ' decreases faster, while the RAAN and Argument of Perigee increase ',...
    ' faster with nrlmsisewhich makes sense given that the nrlmsise model '...
    'is more accurate.'])

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
       % ??????????????????????????????????????????????
    if z > 0
    c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);
    else
    c = 1/2;
    end
 

 end
