function [r_NED,v_NED,v_body]=body2NED(h_sim,v_sim,hea,pitch,bank,par)
%This function transforms point locations from body reference frame (which
%the X axis is directed as the velocity vector, so the aerodynamic frame
%overlaps with the bosy frame, due to the absence of control), to
%local Cartesian coordinates (xNorth,yEast, zDown), given a local coordinate
%system defined by Euler's angles for the reentry phase (heading, pitching,
%banking angles).
%Earth is considered as an ellipsoid
% INPUTs:
%
%v_sim                  [N x 1]                         Velocity vector
%                                                       from simulation
%                                                       done (absolute 
%                                                       magnitude)
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
if nargin ~= 6
    error('Incorrect number of inputs.  See help ecef2enu_12.')
end
if size(h_sim,2) ~= 1
      error('Check the help of this function')
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

r_NED = zeros(3,N);
v_body = zeros(3,N);
v_NED = zeros(3,N);
tic
for j = 1:N
    %Creation of position vector in NED frame
    r_NED(:,j) = [0; 0; -(h_sim(j)+par.Re)];
    %Creation of velocity vector in body frame by doing a proper projection
    v_body(:,j) = [v_sim(j); 0; 0];
    %Creation of cosine matrix from NED2body
    R1 = [cosd(hea(j)), sind(hea(j)), 0; -sind(hea(j)), cosd(hea(j)),0; 0, 0, 1];
    R2 = [cosd(pitch(j)), 0, -sind(pitch(j)); 0, 1, 0; sind(pitch(j)), 0, cosd(pitch(j))];
    R3 = [1, 0, 0; 0, cosd(bank(j)), sind(bank(j)); 0 -sind(bank(j)), cosd(bank(j))];
    %Creation of cosine matrix from body2NED
    R_body2NED = R1*R2*R3;
    v_NED(:,j) = R_body2NED*v_body(:,j);
               
               
end
toc
end
