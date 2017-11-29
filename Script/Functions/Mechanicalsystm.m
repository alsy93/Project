function f = Mechanicalsystm (var)

v = var(1) ;
gamma =var(2);
h = var(3);

f(1,1) = -((rho_s*g)/(2*beta))*v^2 + g*sin(gamma);
f(2,1) = -((rho_s*g)/(2*beta))*(c_l/c_d)*v + (g*cos(gamma)/v)  - v*cos(gamma)/(Re+h);
f(3,1) = -v*sin(gamma);

f = @(t,var) f;
end