function [t, y] = integratorskip(y0,time,par)
     
    function dydt = Mechanicalsystm(t,y,par)

            v = y(1);
            gamma = y(2);
            h = y(3);
            lat = y(4);
            long = y(5);
            hea = y(6);
            
            
            rho = varrho(h);
            
            [Me,Ne] = MeNe(lat,par);

            dydt = zeros(6,1);

            dydt(1,1) = -((rho*par.g)/(2*par.beta))*v^2 - par.g*sind(gamma);
            dydt(2,1) = ((rho*par.g)/(2*par.beta))*par.eff*v + cosd(gamma)*(- (par.g/(v)) + v/(par.Re+ h));
            dydt(3,1) = v*sind(gamma);
            dydt(4,1) = v*cosd(gamma)*cosd(hea)/(Me+ h);
            dydt(5,1) = v*cosd(gamma)*sind(hea)/((Ne+ h)*cosd(lat));
            dydt(6,1) = (v/(par.Re+ h))*cosd(gamma)*sind(hea)*tand(lat);
            
         

    end

% Find the corresponding density at each evaluation


    function rho = varrho(h)
           
        if h > 10 
                   h = h*1e3;
              [rho, ~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
              rho = rho*1e9;                     %UM: [kg/km^3]
        else
            h = h*1e3;
            [~, rho] = atmosisa(h);         %UM: [kg/m^3]
            rho = rho*1e9;                      %UM: [kg/km^3]
        end
    end
   

    function [Me,Ne] = MeNe(lat,par)
        
              Me = (par.Re*(1-par.e^2))./((1-par.e^2.*sin(lat).^2).^(3/2));
              Ne = par.Re./sqrt(1-par.e^2.*sin(lat).^2);
    end

% Try different options

%options = odeset('OutputFcn',@odeplot);
%options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
%options =  odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
%options = odeset('RelTol',1e-8,'AbsTol',1e-9,'OutputFcn',@odeplot,'Stats','on');
% options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
%                @odephas3,'MaxStep',1);
%options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
%                  @odeplot,'OutputSel',[1 2 3 4],'Stats','on','InitialStep',1e-20,...
%                  'Refine',25);

options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
                 @odeplot,'OutputSel',[1 2 3 4 5 6],'Stats','on',...
                 'Events',@gamma_eventskip);
% options15s =  odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
%                  @odeplot,'OutputSel',[1 2 3 4 5 6],'Stats','on','InitialStep',1e-05,...
%                  'MStateDependence','strong','MvPattern','S');
figure()
[t1,y1]= ode45(@Mechanicalsystm,time,y0,options,par);
% figure(3)
% [t3,y3]= ode15s(@Mechanicalsystm,time,y0,options15s,par);
% figure(4)
% [t4,y4]= ode23tb(@Mechanicalsystm,time,y0,options15s,par);
% [t5,y5]= ode23s(@Mechanicalsystm,time,y0,options15s,par);



%Number of function evaluations

% fprintf('No. points = %d, \t fcount = %d \n', size(sol1.y,2), sol1.stats.nfevals) ;
% fprintf('No. points = %d, \t fcount = %d \n', size(sol2.y,2), sol2.stats.nfevals) ;
% fprintf('No. points = %d, \t fcount = %d \n', size(sol3.y,2), sol3.stats.nfevals) ;
%fprintf('No. points = %d, \t fcount = %d \n', size(sol4.y,2), sol4.stats.nfevals) ;

t = t1;y = y1;
% t = t2;y = y2;
% t = t3;y = y3;
%  t = t4;y = y4;
[Me,Ne] = MeNe(y(:,4),par);
rho = varrho (y(:,3));
end