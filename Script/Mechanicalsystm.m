function dydt = Mechanicalsystm(t,y,par)

v = y(1);
gamma = y(2);
h = y(3);
lat = y(4);

dydt = zeros(4,1);

dydt(1,1) = -((par.rho_s*par.g)/(2*par.beta))*v^2 + par.g*sind(gamma);
dydt(2,1) = -((par.rho_s*par.g)/(2*par.beta))*par.eff*v + (par.g*cosd(gamma)/v)  - v*cosd(gamma)/(par.Re+h);
dydt(3,1) = -v*sind(gamma);
dydt(4,1) = v*cosd(gamma)/(par.Re + h);

end