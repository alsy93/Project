% SOYUZ DESCENT MODULE RE-ENTRY (FIRST PHASE)
% Mission MS-03 on June 2, 2017


addpath(genpath('..\Functions'))

% Data of the mechanical part

m = 2.9;       %Mass [kg]
l = 2.1;       %Length [m]
d = 2.2;       %Diameter [m]
Aref = 3.8;    %Reference area [m^2]
eff = 0.26;    %Aerodynamic efficiency
c_l = 0.349;   %Lift coefficient
c_d = 1.341;   %Drag coefficient
h_0 = 99.7;    %Altitude at starting point [km]
v_0 = 7.619;   %Velocity magnitude
gamma_i = -1.5;%Initial flight path angle [�]
lat_0 = 34.35; %Initial latitude [�]
long_0 = 45.56; %Initial longitude [�]

%Define the starting point of simulation: is 3 hours after undocking[y m d h m s]
% and it is referred to Moscow time

t_0 = [2017 06 02 16 47 25];
t_f = [2017 06 02 16 55 39];%Time at which the parachute will be opened [s]

% Data of the termical part

Rn = 2.235; %Nose radius [m]
e = 0.85; %Emissivity coefficient
cd = 1.4; %Conductive coefficient [W/(mK)]

%ECI frame has the Gamma point direction and X-Y plane lies on the
%equatorial plane, while the Z axis passes trought the North Pole.
ECI = eye(3);


