function odei = Mechanicalsystm(t,var,varp)
m = 2.9e3;       %Mass [kg]
Vol = 9;       %Volume [m^3]
rho_s = m/Vol; %Density [kg/m^3]
g = 9.81;   %Gravity acceleration [m/s^2]

v = var(1);
gamma = var(2);
h = var(3);

odei = zeros(3,1);
odei(1) = varp(1) -(((rho_s*g)/(2*beta))*v.^2 + g.*sind(gamma));
odei(2) = varp(2) -(((rho_s*g)/(2*beta))*eff.*v + (g.*cosd(gamma)./v)  - v.*cosd(gamma)./(Re+h));
odei(3) = varp(3) - (-v.*sind(gamma));

end