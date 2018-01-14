function [r_ECEF,v_ECEF]=ned2ecef(r_NED,v_NED,lat,long,h,par)
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
if nargin ~= 6
   error('Incorrect number of inputs.  See help enu2ecef_12.')
end
if size(v_NED,1) ~= 3
   error('Check the help of this function')
end
if size(r_NED,1) ~= 3
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

% The position vector transformation from the geodetic system to the ECEF coordinate
% system is an intermediate step in converting the GPS position measurement to the
% local NED coordinate system. Given a point in the geodetic system defined
% with long, lat, ellipsoidal height the trasformation is the following:

r_geodetic2ECEF = [(Ne+h).*cosd(lat).*cosd(long),...
                   (Ne+h).*cosd(lat).*sind(long),...
                   (Ne.*(1-par.e^2)+h).*sind(lat)];
r_geodetic2ECEF = r_geodetic2ECEF';

               
%Creation of empty matrixes
dr_ECEF = zeros(3,N);
v_ECEF = zeros(3,N);

for j = 1:N
    %Creation of cosine matrix from NED2ECEF
    R_NED2ECEF = [-sind(lat(j))*cosd(long(j)), -sind(long(j)), -cosd(lat(j))*cosd(long(j));...
                  -sind(lat(j))*sind(long(j)), cos(long(j)), - cosd(lat(j))*sind(long(j));...
                  cosd(lat(j)), 0, -sind(lat(j))];
    
    dr_ECEF(:,j) = R_NED2ECEF*r_NED(:,j);
    v_ECEF(:,j) = R_NED2ECEF*v_NED(:,j);
end

r_ECEF = dr_ECEF + r_geodetic2ECEF;

end