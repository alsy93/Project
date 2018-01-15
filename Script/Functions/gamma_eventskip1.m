function [value,isterminal,direction] =gamma_eventskip1(t,y,y0,par)

%This is an event function for the altitude which cannot reach the value of
%11 km

gammamin = 1.5;
value = y(2) - gammamin;
isterminal = 1;
direction = 1;
end