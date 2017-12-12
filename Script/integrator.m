function [t, y] = integrator(y0,time,par)




    function dydt = Mechanicalsystm(t,y,par)

            v = y(1);
            gamma = y(2);
            h = y(3);
            lat = y(4);
            
            rho = varrho(h);

            dydt = zeros(4,1);

            dydt(1,1) = -((rho*par.g)/(2*par.beta))*v^2 - par.g*sind(gamma);
            dydt(2,1) = +((rho*par.g)/(2*par.beta))*par.eff*v - (par.g*cosd(gamma)/(v)) - v*cosd(gamma)/(par.Re + h);
            dydt(3,1) = v*sind(gamma);
            dydt(4,1) = v*cosd(gamma)/(par.Re+ h);

    end
% Find the corresponding density at each evaluation

    function rho = varrho(h)
%             [~, ~, ~, ~, ~, rhof, ~, ~, ~, ~, ~, ~,~ ] = atmo(h);
%             rhof = flipud(rhof);
%             rho = rhof(1);

        if h > 20 %after tropopause
              h = h*1e3;
              [rho,~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
              rho = rho*1e9;                    %UM: [kg/km^3]
        else 
              h = h*1e3;
              [rho,~, ~, ~,~]=tropos(h);        %UM: [kg/m^3]
              rho = rho*1e9;                    %UM: [kg/km^3]
        end
    end



% Try different options

%options = odeset;
%options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
%options =  odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
%options = odeset('RelTol',1e-8,'AbsTol',1e-9,'OutputFcn',@odeplot,'Stats','on');
% options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
%                @odephas3,'MaxStep',1);
%options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
%                  @odeplot,'OutputSel',[1 2 3 4],'Stats','on','InitialStep',1e-20,...
%                  'Refine',25);

options = odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
                 @odeplot,'OutputSel',[1 2 3 4],'Stats','on','InitialStep',1e-5);
options15s =  odeset('RelTol',1e-15,'AbsTol',1e-15,'NormControl','on','OutputFcn',...
                 @odeplot,'OutputSel',[1 2 3 4],'Stats','on','InitialStep',1e-05,'MStateDependence','strong','MvPattern','S');
figure(1)
[t1,y1]= ode45(@Mechanicalsystm,time,y0,options15s,par);
figure(2)
[t2,y2]= ode113(@Mechanicalsystm,time,y0,options,par);
figure(3)
[t3,y3]= ode15s(@Mechanicalsystm,time,y0,options15s,par);
% figure(4)
% [t4,y4]= ode23tb(@Mechanicalsystm,time,y0,options15s,par);
% [t5,y5]= ode23s(@Mechanicalsystm,time,y0,options15s,par);



%Number of function evaluations

% fprintf('No. points = %d, \t fcount = %d \n', size(sol1.y,2), sol1.stats.nfevals) ;
% fprintf('No. points = %d, \t fcount = %d \n', size(sol2.y,2), sol2.stats.nfevals) ;
% fprintf('No. points = %d, \t fcount = %d \n', size(sol3.y,2), sol3.stats.nfevals) ;
%fprintf('No. points = %d, \t fcount = %d \n', size(sol4.y,2), sol4.stats.nfevals) ;

t = t3;y = y3;
% t = t2;y = y2;
% t = t3;y = y3;
% t = t4;y = y4;

% Plotting of the results


    figure(5)
    
      
        ax1 = subplot(2,2,1);
        plot(ax1,t,y(:,1))
        title('Velocity')
        grid on
        
        ax2 = subplot(2,2,2);
        plot(ax2,t,y(:,2))
        title('Flight path angle')
        grid on
        
        ax3 = subplot(2,2,3);
        plot(ax3,t,y(:,3))
        title('Altitude')
        grid on
        
        ax4 = subplot(2,2,4);
        plot(ax4,t,y(:,4))
        title('Latitude')
        grid on
        
    figure(6)
        plot(varrho(y(:,3)),y(:,3))
        title('Vehicle deceleration and atmospheric density')
end
