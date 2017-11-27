function [r_NED,v_NED]=ecef2ned_12(r_ECEF,v_ECEF,lat,long,h)
%This function transforms point locations from geocentric Earth-Centered Earth-Fixed
%(ECEF) coordinates (X, Y, Z) to local Cartesian coordinates (xNorth,
%yEast, zDown), given a local coordinate system defined by the geodetic
%coordinates of its origin settled on the body (lat,long,h). 
%Earth is considered as a perfect sphere
% INPUTs:
%
%r_ECEF                 [3 x N]                         Position Vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%v_ECEF                 [3 x N]                         Velocity vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%lat                    [1 x N]                         Latitude expressed 
%                                                       in degrees starting
%                                                       from Equator [-90,90]
%
%long                   [1 x N]                         Latitude expressed 
%                                                       in degrees starting
%                                                       from Equator [0,360]
%
%h                      [1 x N]                         Altitude expressed in km
%                                                       starting from Earth
%                                                       surface.


%Check input
if nargin ~= 5
    error('Incorrect number of inputs.  See help ecef2enu_12.')
end
if size(r_ECEF,1) ~= 3
   error('Check the help of this function');
end
if size(v_ECEF,1) ~= 3
      error('Check the help of this function')
end

%Checking to see if length of ECEF matrix is the same as the length of the time vector
N = size(r_ECEF,2);
if N ~= length(lat) || N ~= length(long) || N ~= length(h)
    error('Column size of ECEF vector not equal to size of latitude, longitude and altitude vectors. Check inputs.')
end

%Creation of empty matrixes
r_NED = zeros(3,N);
v_NED = zeros(3,N);

for j = 1:N 
    r_NED(:,j) = ecef2ned(r_ECEF(1,j),r_ECEF(2,j),r_ECEF(3,j),lat(j),long(j),h(j),referenceSphere);
    v_NED(:,j) = ecef2ned(v_ECEF(1,j),v_ECEF(2,j),v_ECEF(3,j),lat(j),long(j),h(j),referenceSphere);

end