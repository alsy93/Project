%syms 'v' 'gamma' 'h' 'lat'
% v=var(1);
% gamma=var(2);
% h=var(3);
% 
% dvdt =( -((rho_s*g)/(2*beta))*v.^2 + g.*sind(gamma));
% dgammadt = (-((rho_s*g)/(2*beta))*eff.*v + (g.*cosd(gamma)./v)  - v.*cosd(gamma)./(Re+h))*10^(7);
% dhdt = (-v.*sind(gamma));
% dlatdt = v.*cosd(gamma)./(Re+h)*10^(4);
% 
% Mechanicalsystm=@(t,var)[dvdt;  dgammadt; dhdt; dlatdt];
% options =  odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
% [t,varout] = ode45(Mechanicalsystm,time,var_0,options);


% t_0 = 0;
% var_0 = [v_0 gamma_0 h_0];
% varp_0 = [(v_f-v_0)/delta_t (gamma_f-gamma_0)/delta_t (h_f-h_0)/delta_t];
% 
% v = var(1);
% gamma = var(2);
% h = var(3);
% 
% %odei = zeros(3,1);
% 
% % odei(1) = varp(1) -(((rho_s*g)/(2*beta))*v.^2 + g.*sind(gamma));
% % odei(2) = varp(2) -(((rho_s*g)/(2*beta))*eff.*v + (g.*cosd(gamma)./v)  - v.*cosd(gamma)./(Re+h));
% % odei(3) = varp(3) - (-v.*sind(gamma));
% odei(1) = ((rho_s*g)/(2*beta))*v.^2 + g.*sind(gamma);
% odei(2) =((rho_s*g)/(2*beta))*eff.*v + (g.*cosd(gamma)./v)  - v.*cosd(gamma)./(Re+h);
% odei(3) = -v.*sind(gamma);
% Mechanicalsystm = @(t,var) odei;


% [var_1,varp_1] = decic(Mechanicalsystm,t_0,var_0,1,varp_0,[]);
% options =  odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
% [t,varout] = ode15i(Mechanicalsystm,time,var_1,varp_1,options);
% 
% 
% 