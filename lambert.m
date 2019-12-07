%Lacey Davis
%Lambert's problem for orbits 

function [v1 v2] = lambert(r1,r2, delta_t_given)
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


 
