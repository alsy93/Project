function [v_NED]=body2NED(v_sim,hea,pitch,bank)
%This function transforms point locations from body reference frame to
%local Cartesian coordinates (xNorth,yEast, zDown), given a local coordinate
%system defined by Euler's angles for the reentry phase (heading, pitching,
%banking angles).
%Earth is considered as a perfect sphere
% INPUTs:
%
%v_sim                  [N x 1]                         Velocity vector
%                                                       from simulation
%                                                       done (magnitude)
%
%hea                    [N x 1]                         Heading angle
%                                                       expressed in degree
%                                                       starting from the X
%                                                       body axis in X-Y
%                                                       plane
%
%pitch                  [N x 1]                         Pitch angle
%                                                       expressed in degree
%                                                       starting from the X
%                                                       body axis in X-Z
%                                                       plane
%
%bank                   [N x 1]                         Banking angle
%                                                       expressed in degree
%                                                       starting from Y
%                                                       body axis in Y-Z
%                                                       plane


%Check inputs
if nargin ~= 4
    error('Incorrect number of inputs.  See help ecef2enu_12.')
end
if size(v_sim,2) ~= 1
      error('Check the help of this function')
end
if size(hea,2) ~= 1
      error('Check the help of this function')
end
if size(pitch,2) ~= 1
      error('Check the help of this function')
end
if size(bank,2) ~= 1
      error('Check the help of this function')
end

%Checking to see if length of ECEF matrix is the same as the length of the time vector
N = size(v_sim,2);
if N ~= length(hea) || N ~= length(pitch) || N ~= length(bank)
    error('Column size of velocity body vector not equal to size of heading, pitching and banking vectors. Check inputs.')
end

%

%Creation of empty matrixes
v_NED = zeros(3,N);

for j = 1:N 
    r_NED(:,j) = ecef2ned(r_ECEF(1,j),r_ECEF(2,j),r_ECEF(3,j),lat(j),long(j),h(j),referenceSphere);
    v_NED(:,j) = ecef2ned(v_ECEF(1,j),v_ECEF(2,j),v_ECEF(3,j),lat(j),long(j),h(j),referenceSphere);

end