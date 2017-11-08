function UTC_s = localtime2UTC(time,Timezone)
%This function provides to calculate UTC time starting from local mean time.
% INPUTs:
%
% date:   [year, month, day, hour, minutes, seconds] Local mean time 
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
% Timezone: 1x1       Is the timezone of the place, for example for Moscow is
%                     +3 (UTC/GMT +3 hours).The sign is important!
% OUTPUTs:
%
%UTC_s      1x1             [seconds]  An array providing the UTC in seconds
UTC = datetime(time,'ConvertFrom',Timezone);
UTC_s = datetime(UTC,'ConvertFrom','epochtime','Epoch','2010-3-21-00-00-00');             