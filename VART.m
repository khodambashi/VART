%% Introduction
% This code solve the Variable residence time–based model (VART) for uniform
% channel and flow  parameters with pdepe function in MATLAB.

% This code available at :
%                      https://github.com/khodambashi/VART.git

% For more information : 
%                      https://www.mathworks.com/help/matlab/ref/pdepe.html

% sajad Khodambashi Emami, Ver 1.0, March 2022
%% 
close all
clear
clc

%% Difine Models Parameters
global Ks qL Aadv Tv Q CL A Ch Ds lamda qh
Ch = 0;
Ds = 1 ;




Ks = 0.24; % Fickian dispersion coefficient excluding the transient storage effect[m^2/s]

qL = 0.001 ;  %  the lateral inflow rate on a per length basis (m^3/s-m)

Aadv = 3 ;   % area [L^2] of advection-dominated storage zone.

Tv = 3000 ; %  storage zone exchange coefficient [1/s]

%Adiff = 1 ; % the cross-sectional area of the storage zone [m^2]


CL = 0.5 ; % solute concentration [M/L^3] in lateral inflow.

lamda = 0 ; % dimensionless parameter [].


qh = 2 ; %hyporheic flow gain/loss rate per unit channel length [L2/T].
 %======================= Calculate velocity ==============================
Q = 1 ; % The volumetric flow rate [m^3/s]
A = 12 ; % the main channel cross-sectional area [m^2]

U = Q/A ; % Velocity [m/s]
%% Spatial Mesh
Xmin = 0 ;                % pollutant injection site
L = 281 ;                 % Length of river [m]
NX = 282 ;                % Number of spatial subdivisions
X = linspace(Xmin,L,NX) ; % Subdivision of the Spatial domain.
dX = L/(NX-1) ;           % space-step

%% Temporal Mesh
Tmin = 0 ;                   % initial time
Tmax = (L/U)*1.5  ;          % lengh of simulation [s]
NT = 201 ;                   % Number of temporal subdivisions   
T = linspace(Tmin,Tmax,NT) ; % Subdivision of the temporal domain.
dT = Tmax/(NT-1) ;           % time-step
%% Time Pattern of Pollutant Injection
global t c
c = [10 ,10 ,0 ,0];
t = [0 ,100 ,150 ,10000];
%%  Symmetry constant
m = 0 ;

%% Main Solve Function
C1 = pdepe(m,@VART1func,@VARTic,@VARTbc,X,T);
C2 = pdepe(m,@VART2func,@VARTic,@VARTbc,X,T);

 T(T>=(Tv)) = [];
 ZZ = length(T) +1 ;
 

ZZZ = ZZ + 1 ;

C(1:ZZ,:,1) = C1(1:ZZ,:,1) ;
C(ZZZ:NT,:,1) = C2(ZZZ:NT,:,1) ; 



%% Plot Pollutant Transport in River

for i=1:NT
    
   plot(X,C(i,:,1),'k','LineWidth',1) 
   title('Pllutant Transport in River')
   xlabel('River Route [m]')
   ylabel ('Concentration [ppm]')
   box on
   grid on
   
    axis([0 L 0 max(c)])
    pause(0.0412)
end
close all
%% Difine Equation for pdepe

function [c,f,s] = VART1func(~,~,u,dudx)
global Ks qL Aadv Tv Q CL A Ch lamda qh

c = [1
    1] ;
   

f =  [  Ks * dudx(1)
    0] ;


s = [ (1/Tv)*( (Aadv)/A ) * ( lamda * u(2) - u(1) ) - (Q/A) * dudx(1) + (qL/A) * ( CL - u(1) ) 
   (1/Tv)* ( u(1) - u(2) )+ (qh/(Aadv))*(Ch - u(2)) ] ;

end

function [c,f,s] = VART2func(~,T,u,dudx)
global Ks qL Aadv Q CL A Ch Ds lamda qh

c = [1
    1] ;
   

f =  [  Ks * dudx(1)
    0] ;
Adif = 1;%4*pi*Ds*(T-Tv);

s = [ (1/T+0.000001)*( (Aadv+Adif)/A ) * ( lamda * u(2) - u(1) ) - (Q/A) * dudx(1) + (qL/A) * ( CL - u(1) ) 
   (1/T+0.000001)* ( u(1) - u(2) )+ (qh/(Aadv+Adif))*(Ch - u(2)) ] ;
end

%% Initial Condition
function u0 = VARTic(~)

   u0 = [0
       0];

end

%% Boundray Condition
function [pl,ql,pr,qr] = VARTbc(~,ul,~,~,T)
global t c Ks
C_in = interp1(t,c,T) ;

   % Left Boundary Conditions (time pattern of pollutant injection)
      pl = [ ul(1) - C_in
       ul(2) ];
   
    ql = [0 
        0];
 
   % Right Boundary Conditions ()
    pr = [ 0
        0 ];
    qr = [1/Ks 
        1 ];
 
end

%%-------------------------------------------------------------------------
