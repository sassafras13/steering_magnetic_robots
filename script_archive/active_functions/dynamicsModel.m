% Dynamics Model of Swimmer in Magnetic Gradient Field
% Emma Benjaminson
% 17 Jan 2020

%%%%%%%%
% states
%%%%%%%%
syms x y theta xdot ydot thetadot
X = [x y theta xdot ydot thetadot] ; 

%%%%%%%%%%
% controls
%%%%%%%%%%
syms I1 I2
U = [I1 I2] ; % [A] current in coils 1 and 2

%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%

% microswimmer properties
w = 0.001 * 1000 ; % [mm] width of macroswimmer
L = 0.01 * 1000 ; % [mm] length of macroswimmer
h = 0.005 * 1000 ; % [mm] height of macroswimmer
volume = w * L * h ; % [mm3] volume of macroswimmer
mass = 0.001 * 1000 ; % [g] mass of macroswimmer
mag = 4E3 ; % [A/m] magnetization of link 
drag = 1E-1 ; % [?] drag coefficient 
ct = drag ; % drag in the tangential direction 
cn = 2*ct ; % drag in the normal direction

% environment properties
mu0 = 4*pi*(10^(-7)) ; % [N/A^2] magnetic permeability

% coil dimensions 
a = (0.136 / 2) * 1000 ; % [mm] radius of coils
d = sqrt(3)*(a/2) ; % [mm] Maxwell separation of coils
nturns = 320 ; % number of turns in coils
coils = [I1 ; I2 ; a ; d ; nturns] ; 
Bmax = 0.01 ; % [T] highest value of magnetic field allowed
gBmax = 10 ; % [T/mm] highest value of magnetic field gradient allowed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forces acting on the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force due to magnetic field gradient
[gBx, gBy, gB] = magGradientCoilPoint(coils,mu0,x,y) ; % magnetic field gradient at a point
gBangle = atan2(gBy,gBx) ; 
Fbx = abs(mag*volume*gB*cos(theta - gBangle)) * cos(gBangle) ;
Fby = abs(mag*volume*gB*cos(theta - gBangle)) * sin(gBangle) ; 

% force due to hydrodynamic drag
velAngle = atan2(ydot, xdot) ; 
if (velAngle - gBangle) < pi/4 
    c = -1 ; 
else
    c = 1 ; 
end

that = [cos(theta) ; sin(theta)] ; % tangential vector in body frame
nhat = [-sin(theta) ; cos(theta)] ; % normal vector in body frame

Fht = ct * [xdot ; ydot]' * that ; % tangential force
Fhn = cn * [xdot ; ydot]' * nhat ; % normal force

magFh = sqrt( (Fht^2) + (Fhn^2) ) ; 

Fhx = c * magFh * cos(gBangle) ;
Fhy = c * magFh * sin(gBangle) ;  

% torque due to hydrodynamic drag
p0dot = [xdot ; ydot] ; 
Th = -( ((L^2)/4) * cn * p0dot' * nhat ) - ( ((L^3)/12) * cn * thetadot )   ;

% torque due to misalignment between magnetic field and macroswimmer
% magnetization vector
[Bx, By, B] = magFieldCoilPoint(coils,mu0,x,abs(y),Bmax) ; 
Bangle = atan2(By,Bx) ;
Tb = mag * volume * B * sin(Bangle - theta)  ;

% moment inertia of a rectangle rotating about z axis
I = (1/12) * mass * ((w^2) + (L^2)) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear functions for state variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xddot
xddot = (1/mass) * (Fbx + Fhx) ; 

% yddot
yddot = (1/mass) * (Fby + Fhy) ; 

% thetaddot
thetaddot = (1/I) * (Tb - abs(Th)) ;

XDOT = [xddot ; yddot ; thetaddot] ; 

%%%%%%%%%%%%%%%%%%%%%%%%
% sample operating point
%%%%%%%%%%%%%%%%%%%%%%%%
xstar = [19.3311752229792,25.1726708050330,0,0.162582204699066,-0.0279206444823965,0] ; 
ustar = [-2 ; 0] ; % [A]

%% Function Definitions

function [gBx, gBy, gB] = magGradientCoilPoint(coils,mu0,x,y)   
    % magGradientCoil Calculate the magnetic field gradient created by a
    % pair of wire coils with current flowing through them via the
    % Biot-Savart law. 
    %
    % Inputs: 
    %   coils :- data on the coils, including current in each coil, radius,
    %   spacing between coils and number of turns
    %   coils = [I1 ; I2 ; a ; d ; nturns]
    %
    %   mu0 :- magnetic permeability of free space
    %
    %   x, y :- coordinates of the point at which to calculate the field
    %   gradient
    %
    % Outputs: 
    %   gBx :- x-components of magnetic field gradient at each position in
    %   2D space as defined by positionArray
    %
    %   gBy :- y-components of magnetic field gradient at each position in
    %   2D space as defined by positionArray
    %
    %   gB :- magnitude of magnetic field gradient at each position in 2D
    %   space as defined by positionArray
    
    I1 = coils(1,1) ;
    I2 = coils(2,1) ;
    a = coils(3,1) ; % radius 
    d = coils(4,1) ; 
    nturns = coils(5,1) ;  
    
    x1 = x - d ; 
    [gBx1,gBy1] = dBdi(I1,a,mu0,nturns,x1,y) ; 
    
    x2 = x + d ; 
    [gBx2,gBy2] = dBdi(I2,a,mu0,nturns,x2,y) ; 
    
    gBx = gBx1 + gBx2 ; 
    gBy = gBy1 + gBy2 ; 
    gB = sqrt(gBx.^2 + gBy.^2) ; 
    
    gBx = double(gBx) ; 
    gBy = double(gBy) ; 
    gB = double(gB) ; 


end

function [gBx,gBy] = dBdi(I,a,mu0,n,xarray,yarray)
    % dBdi Calculate the gradient of the magnetic field at every point in
    % xarray, yarray. This function uses an embedded series of symbolic
    % equations to solve the gradient of the magnetic field and then
    % converts the equation to a series of values which it outputs for
    % plotting in a quiver plot or equivalent.
    %
    % Inputs: 
    %   I :- current to coil
    %
    %   a :- radius of coil
    %   
    %   mu0 :- magnetic permeability of free space
    %
    %   n :- number of turns in coil
    %
    %   xarray, yarray :- array of points in space for which we want to
    %   calculate gradient
    %
    % Outputs: 
    %   gBx, gBy :- components of the gradient in the x and y directions
    
    syms x y 'real'
    ksq = (4*a*y) / ((a+y)^2 + x^2) ; 
    [K,E] = ellipke(ksq) ; 
    Bx = ((mu0*I*n) / (2*pi)) * (1 / sqrt((y + a)^2 + x^2)) * (((a^2 - y^2 - x^2) / ((a-y)^2 + x^2))*E + K) ;
    By = ((mu0*I*n) / (2*pi)) * (x / (y * sqrt((y+a)^2 + x^2))) * (((a^2 + y^2 + x^2) / ((a-y)^2 + x^2))*E - K) ;
    B = sqrt(Bx^2 + By^2) ;
    dBdx = simplify(diff(B,x)) ; 
    dBdy = simplify(diff(B,y)) ;
    gBx = subs(dBdx,{x,y},{xarray,yarray}) ; 
    gBy = subs(dBdy,{x,y},{xarray,yarray}) ; 
end

function [Bx, By, B] = magFieldCoilPoint(coils,mu0,x,y,Bmax)
    % magFieldCoil Calculate the magnetic field created by a pair of wire
    % coils with current flowing through them via the Biot-Savart law. 
    %
    % Inputs: 
    %   coils :- data on the coils, including current in each coil, radius,
    %   spacing between coils and number of turns
    %   coils = [I1 ; I2 ; a ; d ; nturns]
    % 
    %   mu0 :- magnetic permeability of free space
    %
    %   positionArray :- array of values defining the extent of the space in
    %   x and y that the function will use to calculate magnetic field values
    %   positionArray = [dx ; xmin ; xmax ; ymin ; ymax]
    %
    %   Bmax :- highest value of magnetic field (in T) allowed 
    %
    % Outputs: 
    %   gx :- meshgrid of x values, to be used for plotting results
    %
    %   gy :- meshgrid of y values, to be used for plotting results
    %
    %   Bx :- x-components of magnetic field at each position in 2D space
    %   as defined by positionArray
    %
    %   By :- y-components of magnetic field at each position in 2D space
    %   as defined by positionArray
    %
    %   B :- magnitude of magnetic field at each position in 2D space as
    %   defined by positionArray
    %
    %   Bxcheck :- analytical solution of magnetic field along center axis
    %   between coils. Only valid for Helmholtz coil configuration! 
    
    I1 = coils(1,1) ; % [A] current
    I2 = coils(2,1) ; % [A] current
    a = coils(3,1) ; % [m] radius 
    d = coils(4,1) ; % [m] separation distance
    nturns = coils(5,1) ; % number of turns
        
    x1 = x - d ;
    ksq1 = ksquared(x1,y,a) ;
    [K1,E1] = ellipke(ksq1) ; 
    Bx1 = BxHH(x1,y,mu0,nturns,I1,a,E1,K1) ;
    By1 = ByHH(x1,y,mu0,nturns,I1,a,E1,K1) ;
        
    % second coil
    x2 = x + d ;         
    ksq2 = ksquared(x2,y,a) ;
    [K2,E2] = ellipke(ksq2) ; 
    Bx2 = BxHH(x2,y,mu0,nturns,I2,a,E2,K2) ;
    By2 = ByHH(x2,y,mu0,nturns,I2,a,E2,K2) ;
        
    % sum 
    Bx = Bx1 + Bx2 ; 
    By = By1 + By2 ; 
      
    Bx(Bx>Bmax) = Bmax ; 
    By(By>Bmax) = Bmax ; 
    
    B = sqrt(Bx.^2 + By.^2) ;
    B(B>Bmax) = Bmax ; 

end

function ksq = ksquared(x,y,a) 
    ksq = (4*a*y) / (((a + y)^2) + (x^2)) ; 
end

function Bx = BxHH(x,y,mu0,nturns,I,a,E,K)

    Bx = ((mu0*nturns*I)/(2*pi)) * (1 / sqrt((y+a)^2 + x^2)) * (((a^2 - y^2 - x^2) / ((a-y)^2 + x^2)) * E + K) ; 
    
end

function By = ByHH(x,y,mu0,nturns,I,a,E,K)
    By = ((mu0 * nturns * I) / (2 * pi)) * (x / (y * sqrt((y + a)^2 + x^2))) * (((a^2 + y^2 + x^2) / ((a - y)^2 + x^2)) * E - K) ; 
end





