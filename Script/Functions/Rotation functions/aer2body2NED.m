function [v_NED,v_body]=aer2body2NED(v_sim,hea,pitch,bank)
%This function transforms point locations from body reference frame to
%local Cartesian coordinates (xNorth,yEast, zDown), given a local coordinate
%system defined by Euler's angles for the reentry phase (heading, pitching,
%banking angles).
%Earth is considered as an ellipsoid
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
% OUTPUTs:
%
%v_body                 [3 x N]                         Velocity vector in
%                                                       body reference
%                                                       frame
%
%v_NED                  [3 x N]                         Velocity vector in
%                                                       Nord, East, Down
%                                                       frame


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
N = size(v_sim,1);
if N ~= length(hea) || N ~= length(pitch) || N ~= length(bank)
    error('Column size of velocity body vector not equal to size of heading, pitching and banking vectors. Check inputs.')
end

%Trasformation from aerodynamic to body frame [3 x N] (case of no banking)
%and rotation in the NED frame

v_body = zeros(3,N);
v_NED = zeros(3,N);
tic
for j = 1:N
    %Creation of velocity vector in body frame by doing a proper projection
    v_body(:,j) = [v_sim(j)*cosd(hea(j)) + v_sim(j)*cosd(pitch(j));...
                   v_sim(j)*sind(hea(j)); v_sim(j)*sind(pitch(j))];
    %Creation of cosine matrix from NED2body
    R1 = [cosd(hea(j)), sind(hea(j)), 0; -sind(hea(j)), cosd(hea(j)),0; 0, 0, 1];
    R2 = [cosd(pitch(j)), 0, -sind(pitch(j)); 0, 1, 0; sind(pitch(j)), 0, cosd(pitch(j))];
    R3 = [1, 0, 0; 0, cosd(bank(j)), sind(bank(j)); 0 -sind(bank(j)), cosd(bank(j))];
    
    R_NED2body = R1*R2*R3;
    R_body2NED = R_NED2body';
    
    v_NED(:,j) = R_body2NED*v_body(:,j);
               
               
end
toc
end
