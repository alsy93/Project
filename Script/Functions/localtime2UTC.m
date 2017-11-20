function UTC_s = localtime2UTC(time,timezone)

%This function provides to calculate UTC time starting from local mean time.
% INPUTs:
%
% time:   [year, month, day, hour, minutes, seconds] Local mean time 
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
% Timezone: 1x1      An ISO 8601 character  of the form +HH,mm or -HH,mm;
%                    for example, +3.00, to specify a time zone that is a 
%                    fixed offset from UTC (like Moscow).
% OUTPUTs:
%
%UTC_s      1x1             [seconds]  An array providing the UTC in seconds
%                            elapsed since 00:00:00 1-Jan-1970 UTC as the
%                            Posix definition

local_time = datetime(time);
UTC = local_time - sign(timezone)*hours(timezone);
UTC_s = posixtime(UTC);
end



