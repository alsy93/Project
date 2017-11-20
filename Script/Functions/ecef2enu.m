function [r_ENU,v_ENU]=ecef2enu(r_ECEF,v_ECEF,lat,long,h)
% INPUTs:
%
%r_ECEF                 [3 x N]                         Position Vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%v_ECEF                 [3 x N]                         Velocity vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%to be completed

%Check input
if nargin ~= 5
    error('Incorrect number of inputs.  See help eci2ecef.')
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
    error('Size of ECEF vector not equal to size of latitude, longitude and altitude vectors. Check inputs.')
end

%Creation of empty matrixes
r_ENU = zeros(3,N);
v_ENU = zeros(3,N);

for j = 1:N 
    r_ENU(:,j) = ecef2enu(r_ECEF(1,j),r_ECEF(2,j),r_ECEF(3,j),lat(j),lon(j),h(j),referenceSphere);

end