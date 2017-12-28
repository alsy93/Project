function [value,isterminal,direction] =h_event(t,y,y0,par)

%This is an event function for the altitude which cannot reach the value of
%11 km

hmin = 10.8;
value = y(3)-hmin;
isterminal = 1;
direction = -1;

end