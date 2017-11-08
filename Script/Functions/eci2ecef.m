function [r_ECEF,v_ECEF]=eci2ecef(time,r_ECI,v_ECI)
% La function è per definire la terna fissa con la terra, attraverso la
% rotazione della terna centrata inerziale della terra (ECI)
%
% INPUTs:
% 
% time:   [year, month, day, hour, minutes, seconds] Universal Coordinated Time (UTC) 
%           
%            year: enter a double value that is a whole number greater than
%                    1, such as 2013 
%            month: enter a double value that is a whole number greater than
%                    0, within the range 1 to 12.
%            day:  enter a double value that is a whole number greater than
%                    0, within the range 1 to 31.
%            hour: enter a double value that is a whole number greater than
%                    0, within the range 1 to 24.
%            minutes: enter a double value that is a whole number greater than 
%                     0, within the range 1 to 60.
%            seconds:  enter a double value that is a whole number greater 
%                      than 0, within the range 1 to 60.
%r_ECI                  [3 x N]                         Position Vector
%                                                       [km]       
%                                                       in ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity Vector
%                                                       [km/s]   in
%                                                       ECI coordinate
%                                                       frame of reference
%
%a_ECI                  [3 x N]                         Acceleration Vector
%                                                       [km/s2]   
%                                                       in ECI coordinate
%                                                       frame of reference
%
% OUTPUTs
%
%r_ECEF                 [3 x N]                         Position Vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%v_ECEF                 [3 x N]                         Velocity vector in
%                                                       ECEF coordinate
%                                                       frame of reference
%
%a_ECEF                 [3 x N]                         Acceleration Vector
%                                                       in ECEF coordinate
%                                                       frame of reference
%

%Define the reduction actually most used in astrodynamic: it refers to
%coordinates for January 12, 2000 at 4 hours, 52 minutes, 12.4 seconds and 
%January 12, 2000 at 4 hours, 52 minutes, and 13 seconds.

reduction = 'IAU-2000/2006';

%Define DCIP from Earth Rotation and Reference Systems Service (IERS)
IERS=importdata('IERS.txt')

%Define the director cosine matrix 

dcm = dcmeci2ecef(reduction,time);

%Trasformation in ECEF reference frame

r_ECEF = dcm*[r_ECI(1); r_ECI(2); r_ECI(3)];
v_ECEF = diff(dcm)*[v_ECI(1); v_ECI(2); v_ECI(3)];




end

    
    
