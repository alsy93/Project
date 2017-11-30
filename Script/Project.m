clear all; close all
%dbstop if naninf
%% SOYUZ DESCENT MODULE RE-ENTRY (FIRST PHASE)
% Mission MS-03 on June 2, 2017


addpath(genpath('Functions'))

% Data of the mechanical part
%Body data
m = 2.9e3;       %Mass [kg]
l = 2.1;       %Length [m]
d = 2.2;       %Diameter [m]
Vol = 9;       %Volume [m^3]
rho_s = m/Vol; %Density [kg/m^3]
A_ref = 3.8;    %Reference area [m^2]
%Aerodynamic data
eff = 0.26;    %Aerodynamic efficiency
c_l = 0.349;   %Lift coefficient
c_d = 1.341;   %Drag coefficient
%Simulation data
Re = 6738e3;
g = 9.81;   %Gravity acceleration [m/s^2]      
gamma_0 = -1.5;%Initial flight path angle [�]
gamma_0rad = gamma_0*(pi/180);%Initial flight path angle[rad]
h_0 = 99.7e3;    %Altitude at starting point [m]
h_f = 10.8e3;    %Altitude in which modelling finish and parachute will opened [m]
v_0 = 7.619e3;   %Initial velocity magnitude [m/s]
v_f = 0.216e3;   %Final velocity magnitude [m/s]
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

Rn = 2.235; %Nose radius [m]
e = 0.85; %Emissivity coefficient
cd = 1.4; %Conductive coefficient [W/(mK)]

%% Start modeling and simulation

%ECI frame has the Gamma point direction and X-Y plane lies on the
%equatorial plane, while the Z axis passes trought the North Pole.
ECI = eye(3);

%Create a time discretization
UTC_ref = posixtime(datetime(t_ref));
UTC_t_0 = localtime2UTC(t_0,timezone) - UTC_ref;
UTC_t_f = localtime2UTC(t_f,timezone) - UTC_ref;
delta_t = UTC_t_f - UTC_t_0;

N = 5000;
time = linspace(0,delta_t,N); %almost every second
%time = linspace(UTC_t_0,UTC_t_f,N);

%Build system of ODE's
var_0 = [v_0 gamma_0 h_0];
% syms v, gamma,h,lat
v=var(1);
gamma=var(2);
h=var(3);

dvdt = -((rho_s*g)/(2*beta))*v.^2 + g.*sind(gamma);
dgammadt = -((rho_s*g)/(2*beta))*eff.*v + (g.*cosd(gamma)./v)  - v.*cosd(gamma)./(Re+h);
dhdt = (-v.*sind(gamma));

Mechanicalsystm=@(t,var)[dvdt;  dgammadt; dhdt];
options =  odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
[t,varout] = ode15s(Mechanicalsystm,time,var_0,options);

