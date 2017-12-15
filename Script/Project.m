clc;clear all; close all; hold on
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
g = 9.81e-3;   %Gravity acceleration [km/s^2]      
gamma_0 = -1.5;%Initial flight path angle [�]
gamma_0rad = gamma_0*(pi/180);%Initial flight path angle[rad]

h_0 = 99.7;    %Altitude at starting point [m]
h_f = 10.8;    %Altitude in which modelling finish and parachute will opened [m]
v_0 = 7.619;   %Initial velocity magnitude [m/s]
v_f = 0.216;   %Final velocity magnitude [m/s]
beta = (m*g/(c_d*A_ref)); %Ballistic coefficient ratio [Pa]
%---------------------------------------------------------------------------------
lat_0 = 34.35; %Initial latitude [�]
lat_f = 47.22; %Final latitude [�]
long_0 = 45.56;%Initial longitude [�]
long_f = 69.36;%Final longitude [�]
hea_0 = 47.03; %Initial heading angle [�]
hea_f =  62.80;%Final heading angle [�]

% Define the starting point of simulation: is 3 hours after undocking[y m d h m s]
% and it is referred to Moscow time
% Delta_s = 494 s

t_0 = [2017 06 02 16 47 25];   %Initial Time   
t_f = [2017 06 02 16 55 39];   %Time at which the parachute will be opened [s]
t_ref = [2010 01 01 00 00 00]; %Referring time
timezone = +3.00;              %Moscow timezone

% Data of the termical part

Rn = 2.235;             %Nose radius of stagnation point [m]
emiss = 0.85;           %Emissivity coefficient
cd = 1.4;               %Conductive coefficient [W/(mK)]
sigma = 5.670373e-8;    %Stefan-Boltzmann constant [W/(m^2*K^4)]
C = 1;                  % First constant of Tauber-Sutton
a = 1;                  % Second constant of Tauber-Sutton
b = 1.2;                % Third constant of Tauber-Sutton
d = 8.5;                % Fourth constant of Tauber-Sutton


%======== Start modeling and simulation ================================

%Create a time discretization
UTC_ref = posixtime(datetime(t_ref));
UTC_t_0 = localtime2UTC(t_0,timezone) - UTC_ref;
UTC_t_f = localtime2UTC(t_f,timezone) - UTC_ref;
delta_t = UTC_t_f - UTC_t_0;

% N = 500;
% time = linspace(0,delta_t,N); %almost every second
%time = linspace(UTC_t_0,UTC_t_f,N);
time = [0 delta_t];

% Build system of ODE's

y0.v0 = v_0;
y0.gamma0 = gamma_0;
y0.h0 = h_0;
y0.lat0 = lat_0;
y0.long0 = long_0;
y0.hea0 = hea_0;
y0 = [y0.v0 y0.gamma0 y0.h0 y0.lat0 y0.long0 y0.hea0];

par.rho_s = rho_s;
par.beta = beta;
par.g = g;
par.eff = c_l/c_d;
par.Re = Re;

% Solve mechanical part

[t, y] = integrator(y0,time,par);

%Convert km in m

v = y(:,1).*1e3;
h = y(:,3).*1e3;

% Solve thermal part

parT.Rn = Rn;
parT.emiss = emiss;
parT.sigma = sigma;
parT.cd = cd;
parT.Rn = Rn;
parT.a = a;
parT.b = b;
parT.C = C;
parT.d = d;


%ECI frame has the Gamma point direction and X-Y plane lies on the
%equatorial plane, while the Z axis passes trought the North Pole.
ECI = eye(3);
