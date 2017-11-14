function [r_ECEF,v_ECEF]=eci2ecef(r_ECI,v_ECI,UTC_time)
% La function è per definire la Earth centred fixed frame (ECEF), attraverso la
% rotazione della terna centrata inerziale della terra (ECI)
%
% INPUTs:
%
%UTC_time                   [1 x N]                        A vector providing
%                                                      the UTC in seconds
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
% OUTPUTs:
%
%r_ECEF                 [3 x N]                         Position Vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%v_ECEF                 [3 x N]                         Velocity vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%Check inputs
if nargin ~= 3
    error('Incorrect number of inputs.  See help eci2ecef.')
end
if size(r_ECI,1) ~= 3
   error('Check the help of this function');
end
if size(v_ECI,1) ~= 3
      error('Check the help of this function')
end

%Checking to see if length of ECI matrix is the same as the length of the time vector
N = size(r_ECI,2);
if N ~= length(UTC_time)
    error('Size of ECI vector not equal to size of time vector.  Check inputs.')
end

r_ECEF = zeros(3,N);
v_ECEF = zeros(3,N);

%Definition of Greenwich mean time at t0 (fixed at 21 March 2010 at 00:00:00)
%and of Greenwich hour angle
we = (2*pi/365.26)/(24*3600); %angular velocity of the Earth
UTC_t0 =  posixtime(datetime([2010 03 21 00 00 00]))*ones(1,N);
GHA = zeros(1,N);

%Definition of Greenwich hour angle (Greenwich sideral angle at t0)
GST = GHA + we*(UTC_time-UTC_t0);

%Creation of direction cosine matrix for ECEF
dcm_ECEF2ECI = [cos(GST) -sin(GST) 0;sin(GST) cos(GST) 0; 0 0 1];
dcm_ECI2ECEF = dcm_ECEF2ECI';

%Creation the derivated matrix to find velocity 
dcm_dot_ECEF2ECI = [-sin(GST).*we, - cos(GST).*we, 0; cos(GST).*we, -sin(GST).*we, 0;0, 0, 1];
dcm_dot_ECI2ECEF = dcm_dot_ECEF2ECI';

for j = 1:N  %Iterating thru the number of positions provided by user
             % Rotating the ECI vector into the ECEF frame via the GST angle about the Z-axis
    r_ECEF(:,j) = dcm_ECI2ECEF(GST(j))*r_ECI(:,j);
    v_ECEF(:,j) = dcm_dot_ECI2ECEF(GST(j))*r_ECI(:,j) + dcm_ECI2ECEF(GST(j))*v_ECI(:,j);
end

end
