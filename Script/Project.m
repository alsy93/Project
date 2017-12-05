clear all; close all
%dbstop if naninf
%% SOYUZ DESCENT MODULE RE-ENTRY (FIRST PHASE)
% Mission MS-03 on June 2, 2017


addpath(genpath('Functions'))

% Data of the mechanical part
%Body data
m = 2.9e3;       %Mass [kg]
l = 2.1e-3;       %Length [km]
d = 2.2e-3;       %Diameter [km]
Vol = 9e-9;       %Volume [km^3]
rho_s = m/Vol; %Density [kg/km^3]
A_ref = 3.8e-6;    %Reference area [km^2]
%Aerodynamic data
eff = 0.26;    %Aerodynamic efficiency
c_l = 0.349;   %Lift coefficient
c_d = 1.341;   %Drag coefficient
%Simulation data
Re = 6738;     %Radius of the Earth [km]
g = 9.81e-3;   %Gravity acceleration [m/s^2]      
gamma_0 = -1.5;%Initial flight path angle [�]
gamma_0rad = gamma_0*(pi/180);%Initial flight path angle[rad]

h_0 = 99.7;    %Altitude at starting point [m]
h_f = 10.8;    %Altitude in which modelling finish and parachute will opened [m]
v_0 = 7.619;   %Initial velocity magnitude [m/s]
v_f = 0.216;   %Final velocity magnitude [m/s]
beta = (m*g/(c_d*A_ref)); %Ballistic coefficient ratio [kPa*m/s^2]
%---------------------------------------------------------------------------------
lat_0 = 34.35; %Initial latitude [�]
lat_0rad = lat_0*(pi/180);
lat_f = 47.22; %Final latitude [�]
long_0 = 45.56;%Initial longitude [�]
long_f = 69.36;

% Define the starting point of simulation: is 3 hours after undocking[y m d h m s]
% and it is referred to Moscow time
% Delta_s = 494 s

t_0 = [2017 06 02 16 47 25];   %Initial Time   
t_f = [2017 06 02 16 55 39];   %Time at which the parachute will be opened [s]
t_ref = [2010 01 01 00 00 00]; %Referring time
timezone = +3.00;              %Moscow timezone

% Data of the termical part

Rn = 2.235e-3; %Nose radius [m]
e = 0.85; %Emissivity coefficient
cd = 1.4; %Conductive coefficient [W/(mK)]

%======== Start modeling and simulation ================================

%Create a time discretization
UTC_ref = posixtime(datetime(t_ref));
UTC_t_0 = localtime2UTC(t_0,timezone) - UTC_ref;
UTC_t_f = localtime2UTC(t_f,timezone) - UTC_ref;
delta_t = UTC_t_f - UTC_t_0;

N = 5000;
time = linspace(0,delta_t,N); %almost every second
%time = linspace(UTC_t_0,UTC_t_f,N);

% %Build system of ODE's

y0.v0 = v_0;
y0.gamma0 = gamma_0;
y0.h0 = h_0;
y0.lat0 = lat_0;
y0 = [y0.v0 y0.gamma0 y0.h0 y0.lat0];

par.rho_s = rho_s;
par.beta = beta;
par.g = g;
par.eff = eff;
par.Re = Re;

% Try different options

%options = odeset;
%options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
%options =  odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
%options = odeset('RelTol',1e-8,'AbsTol',1e-9,'OutputFcn',@odeplot,'Stats','on');
%options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
                 %@odephas3,'MaxStep',1e-2*abs(delta_t));
% options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
%                  @odeplot,'OutputSel',[1 2 3 4],'Stats','on','InitialStep',1e-20,...
%                  'Refine',25);

options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
                 @odeplot,'OutputSel',[1 2 3 4],'Stats','on','InitialStep',1e-20,...
                 'Events',@);

[t,y] = ode15s(@Mechanicalsystm,time,y0,options,par);



%ECI frame has the Gamma point direction and X-Y plane lies on the
%equatorial plane, while the Z axis passes trought the North Pole.
ECI = eye(3);
