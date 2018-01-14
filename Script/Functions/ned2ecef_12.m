function [v_ECEF]=ned2ecef(v_NED,lat,long,h,par)
%This function transforms point locations from local Cartesian coordinates (xNorth,
% yEast, zDown) to geocentric Earth-Centered Earth-Fixed
%(ECEF) coordinates (X, Y, Z), given a local coordinate system defined by the geodetic
%coordinates of its origin settled on the body (lat,long,h)
%Earth is considered as an ellipsoid with the following properties:
%            LengthUnit: 'kilometer'
%         SemimajorAxis: 6378.137
%         SemiminorAxis: 6356.75231414036
%     InverseFlattening: 298.257222101
%          Eccentricity: 0.0818191910428158
% INPUTs:
%
%v_ECEF                 [3 x N]                         Velocity vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%lat                    [N x 1]                         Latitude expressed 
%                                                       in degrees starting
%                                                       from Equator [-90,90]
%
%long                   [N x 1]                         Latitude expressed 
%                                                       in degrees starting
%                                                       from Equator [0,360]
%
%h                      [N x 1]                         Altitude expressed in km
%                                                       starting from Earth
%                                                       surface.


%Check input
if nargin ~= 5
   error('Incorrect number of inputs.  See help enu2ecef_12.')
end
if size(v_NED,1) ~= 3
   error('Check the help of this function')
end

%Checking to see if length of ECEF matrix is the same as the length of the time vector
N = size(v_NED,2);
if N ~= length(lat) || N ~= length(long) || N ~= length(h)
   error('Column size of ENU vector not equal to size of latitude, longitude and altitude vectors. Check inputs.')
end

    function [Me,Ne] = MeNe(lat,par)
        
              Me = (par.Re*(1-par.e^2))./((1-par.e^2.*sin(lat).^2).^(3/2));
              Ne = par.Re./sqrt(1-par.e^2.*sin(lat).^2);
    end

[~,Ne] = MeNe(lat,par);

%Creation of empty matrixes
v_ECEF = zeros(3,N);


for j = 1:N 
    v_ECEF(:,j) = ned2ecef(v_NED(1,j),v_NED(2,j),v_NED(3,j),lat(j),long(j),h(j),grs80);
end
end