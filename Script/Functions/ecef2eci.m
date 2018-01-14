function [r_ECI,v_ECI]=ecef2eci(r_ECEF,v_ECEF,UTC_time)
% La function è per definire la Earth centred fixed frame (ECEF), attraverso la
% rotazione della terna centrata inerziale della terra (ECI)
%
%INPUTs:
%
%UTC_time               [1 x N]                        A vector providing
%                                                      the UTC in seconds
%
%r_ECEF                 [3 x N]                         Position Vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%v_ECEF                 [3 x N]                         Velocity vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%OUTPUTs:
%
%r_ECI                  [3 x N]                         Position Vector
%                                                       [km]       
%                                                       in ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity Vector
%                                                       [km/s]   in
%                                                       ECI coordinate
%                                                       frame of reference

%Check input
if nargin ~= 3
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
if N ~= length(UTC_time)
    error('Size of ECEF vector not equal to size of time vector.  Check inputs.')
end

%Definition of Greenwich mean time at t0 (fixed at 2 June 2017 at 13:47:25)
%and of Greenwich hour angle
we = (2*pi/365.26)/(24*3600); %angular velocity of the Earth
UTC_t0 =  posixtime(datetime([2010 03 21 13 47 25]))*ones(1,N);
GHA = zeros(1,N);

%Definition of Greenwich hour angle (Greenwich sideral angle at t0)
GST = GHA + we*(UTC_time-UTC_t0);

%Creation of empty matrixes
r_ECI = zeros(3,N);
v_ECI = zeros(3,N);

for j = 1:N  %Iterating thru the number of positions provided by user
             % Rotating the ECI vector into the ECEF frame via the GST angle about the Z-axis
    
    %Creation of direction cosine matrix for ECEF
    dcm_ECEF2ECI = [cos(GST(j)) -sin(GST(j)) 0;...
                    sin(GST(j)) cos(GST(j)) 0;...
                    0 0 1];

    %Creation the derivated matrix to find velocity 
    dcm_dot_ECEF2ECI = [-sin(GST(j))*we, - cos(GST(j))*we, 0;...
                        cos(GST(j))*we, -sin(GST(j))*we, 0; 0, 0, 1];
         
    r_ECI(:,j) = dcm_ECEF2ECI*r_ECEF(:,j);
    v_ECI(:,j) = dcm_dot_ECEF2ECI*r_ECEF(:,j) + dcm_ECEF2ECI*v_ECEF(:,j);
end
end
