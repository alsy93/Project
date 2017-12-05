function dydt = Mechanicalsystm(t,y,par)
% m = 2.9e3;       %Mass [kg]
% Vol = 9;       %Volume [m^3]
% rho_s = m/Vol; %Density [kg/m^3]
% g = 9.81;   %Gravity acceleration [m/s^2]

v = y(1);
gamma = y(2);
h = y(3);
lat = y(4);

dydt = zeros(4,1);

dydt(1,1) = -((par.rho_s*par.g)/(2*par.beta))*v^2 + par.g*sind(gamma);
dydt(2,1) = -((par.rho_s*g)/(2*par.beta))*par.eff*v + (par.g*cosd(gamma)/v)  - v*cosd(gamma)/(par.Re+h);
dydt(3,1) = -v*sind(gamma);
dydt(4,1) = v*cosd(gamma)/(par.Re + h);

end