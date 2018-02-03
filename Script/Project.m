clc;
clear all;
close all;

%==========================================================================
% SOYUZ DESCENT MODULE RE-ENTRY (FIRST PHASE)
% Group 12 : Aloisia Russo and Giulia Bortolato
% Mission MS-05 on December 14th, 2017
% Matlab version R2017b
%==========================================================================

% Options settings

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultTextFontSize',12);
set(0,'DefaultLegendFontSize',5);
set(0,'defaultAxesFontSize',12);
addpath(genpath('Functions'))

% Data of the mechanical part

%Body data
m = 2.9e3;          %Mass [kg]
l = 2.1e-3;         %Length [km]
d = 2.2e-3;         %Diameter [km]
Vol = 4;            %Volume [m^3]
A_ref = 3.8e-6;     %Reference area [km^2]

%Aerodynamic data
c_l = 0.349;   %Lift coefficient
c_d = 1.341;   %Drag coefficient
eff = c_l/c_d; %Aerodynamic efficiency

%Simulation data
Re = 6738;              %Radius of the Earth at equator [km]
Rp = 6357;              %Radius of the Earth at the poles [km]
e = sqrt(Re^2-Rp^2)/Re; %Eccentricity due to Earth oblateness
g = 9.81e-3;            %Gravity acceleration [km/s^2] 
om = 2*pi*Re*1e3/(24*3600); %Velocity of Earth rotation [rad/s]

h_0 =99.5;              %Altitude at starting point [km]
h_f = 10.7;             %Altitude in which modelling finish and parachute will opened [km]
v_0 = 7.618;            %Initial velocity magnitude [km/s]
v_f = 0.218;            %Final velocity magnitude [km/s]
beta = (m/(c_d*A_ref)); %Ballistic coefficient ratio [kg/(km^2)]
alpha = (m/(c_l*A_ref));%Glide coefficient ratio [kg/(km^2)]

%Angles data
gamma_0 = -1.5;%Initial flight path angle [�]
lat_0 = 34.06; %Initial latitude [�]
lat_f = 47.19; %Final latitude [�]
long_0 = 45.26;%Initial longitude [�]
long_f = 69.34;%Final longitude [�]
hea_0 = 45;    %Initial heading angle [�]
hea_f = 62.80; %Final heading angle [�]

% Define the starting point of simulation: is 3 hours after undocking[y m d h m s]
% and it is referred to Moscow time
% Delta_s = 494 s

t_0 = [2017 12 14 11 15 07];   %Initial Time   
t_f = [2017 12 14 11 23 28];   %Time at which the parachute will be opened [s]
t_ref = [2010 01 01 00 00 00]; %Referring time
timezone = +3.00;              %Moscow timezone
%---------------------------------------------------------------------------------
%Data of the termical part

A_refT = 3.8;           %Reference area [m^2]
Rn = 2.235e2;           %Nose radius of stagnation point [cm]
emiss = 0.85;           %Emissivity coefficient
cd = 1.4;               %Conductive coefficient [W/(mK)]
sigma = 5.670373e-8;    %Stefan-Boltzmann constant [W/(m^2*K^4)]
C = 1;                  %First constant of Tauber-Sutton
a = 1;                  %Second constant of Tauber-Sutton
b = 1.6;                %Third constant of Tauber-Sutton
d = 8.5;                %Fourth constant of Tauber-Sutton
Ks = 1.7415e-4;         %Costant of Sutton-Graves for Earth
 
%======== Start modeling and simulation ================================

%Create a time discretization
UTC_ref = posixtime(datetime(t_ref));
UTC_t_0 = localtime2UTC(t_0,timezone) - UTC_ref;
UTC_t_f = localtime2UTC(t_f,timezone) - UTC_ref;
delta_t = UTC_t_f - UTC_t_0;

time = [0 delta_t];

% Build system of ODE's

y0.v0 = v_0;
y0.gamma0 = gamma_0;
y0.h0 = h_0;
y0.lat0 = lat_0;
y0.long0 = long_0;
y0.hea0 = hea_0;
y0 = [y0.v0 y0.gamma0 y0.h0 y0.lat0 y0.long0 y0.hea0];


par.beta = beta;
par.alpha = alpha;
par.g = g;
par.Re = Re;
par.e = e;
par.om = om;
par.cl = c_l;
par.cd = c_d;
par.S = A_ref;
par.m = m;

% Solve mechanical part


[t, y, bank, TimeBank] = integrator_MOD(y0,time,par);


%Convert km in m

v = y(:,1);
gamma = y(:,2);             %flight path angle which corresponds to the pitch angle
h = y(:,3);
h_m = h.*1e3;
lat = y(:,4);
long = y(:,5);
hea = y(:,6);               %heading angle



%% Solve thermal part
parT.Rn = Rn;
parT.emiss = emiss;
parT.sigma = sigma;
parT.cd = cd;
parT.Rn = Rn;
parT.a = a;
parT.b = b;
parT.C = C;
parT.d = d;
parT.Ks = Ks;
parT.A_ref = A_refT;
parT.Vol = Vol;

q_conv = convective_flux(v,h_m,parT);
[T_w, q_rad_TS] = wall_temperature(v,h_m,q_conv,t,parT);

% Change dimension of convective and radiative flux

q_rad_TS_m = q_rad_TS.*1e4;   %[W/m^2]
q_conv_m = q_conv.*1e4;       %[W/m^2]

[T_cab,x] = thermal_shield(h,v,T_w,q_rad_TS_m,q_conv_m,parT,t);

%% Flight envelope

flight_envelope(h,v,T_w,par,parT);
%% Ground track

ground_track(long,lat);


%% Integrator analysis

% Consider a default time step

% For ODE113
[t113_0, y113_0, bank113_0, timeBank113_0] = integrator_MOD113_0(y0,time,par);
plot_integrator_0(t,y,bank,TimeBank,t113_0,y113_0,bank113_0,timeBank113_0);
    

%By using a Max step of 0.3

[t45_1, y45_1, bank45_1,timeBank45_1] = integrator_MOD45_1(y0,time,par);
[t113_1, y113_1,bank113_1,timebank113_1] = integrator_MOD113_1(y0,time,par);

Vargout = plot_integrator_1(t45_1, y45_1, bank45_1,timeBank45_1,t113_1, y113_1,bank113_1,timebank113_1);


%% Sensitivity analysis
% By changing flight path angle

%1) gamma = -0.75;

gamma0075 = -0.75;
y0_075 = [v_0 gamma0075 h_0 lat_0 long_0 hea_0];
[t_075, y_075, bank_075, TimeBank_075] = integrator_MOD(y0_075,time,par);

v_075 = y_075(:,1);
gamma_075 = y_075(:,2);             %flight path angle which corresponds to the pitch angle
h_075 = y_075(:,3);
h_m_075 = h_075.*1e3;
lat_075 = y_075(:,4);
long_075 = y_075(:,5);
hea_075 = y_075(:,6); 

q_conv_075 = convective_flux(v_075,h_m_075,parT);
[T_w_075, q_rad_TS_075] = wall_temperature(v_075,h_m_075,q_conv_075,t_075,parT);

ground_track(long_075,lat_075);

%1) gamma = -3.5;

gamma035 = -3.5;
y0_35 = [v_0 gamma035 h_0 lat_0 long_0 hea_0];
[t_35, y_35, bank_35, TimeBank_35] = integrator_MOD(y0_35,time,par);

v_35 = y_35(:,1);
gamma_35 = y_35(:,2);             %flight path angle which corresponds to the pitch angle
h_35 = y_35(:,3);
h_m_35 = h_35.*1e3;
lat_35 = y_35(:,4);
long_35 = y_35(:,5);
hea_35 = y_35(:,6); 

q_conv_35 = convective_flux(v_35,h_m_35,parT);
[T_w_35, q_rad_TS_35] = wall_temperature(v_35,h_m_35,q_conv_35,t_35,parT);

ground_track(long_35,lat_35);




flight_envelope_sensitivity(h_075,v_075,T_w_075,h_35,v_35,T_w_35,par,parT);





