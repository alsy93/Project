function [T_w, q_rad_TS] = wall_temperature(v,h,q_conv,t,parT)

% Function to calculate the wall temperature at the stagnation point.
% It comes from an equality of two equations, the analytic Stefan-Boltzman
% equation for radiative heat flux and the sum of statistical Sutton-Tauber
% one and the convective heat flux.

% Find corresponding density 
    function rho = varrho(h)
           
        if h > 10000
           [rho, ~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
           
        else
           [~, ~, ~, rho] = atmosisa(h);      %UM: [kg/m^3]
           
        end
    end
rho = varrho(h);

% Start finding corresponding wall temperature 
q_rad_TS = (((parT.Rn*1e-2)^parT.a) .* (rho.^parT.b) .*(v.^(parT.d))); %[W/cm^2]
q_rerad = (q_conv + q_rad_TS).*1e4;                            
T_w = (q_rerad ./ (parT.sigma * parT.emiss)).^(1/4);

% Plotting  of results

figure()
        ax1 = subplot(2,3,1);legend('Location','NW'); hold on; 
        plot(ax1,v,q_rad_TS,'-b','LineWidth',2,'DisplayName','Radiative heat transfer')
        plot(ax1,v,q_conv,'-r','LineWidth',2,'DisplayName','Convectie heat transfer')
%         ax = gca;
%         ax.YScale = 'log';
        xlabel('Velocity $[\frac{km}{s}]$')
        ylabel('Heat transfer $[\frac{W}{cm^2}]$')
        grid on;grid minor
        
        ax2 = subplot(2,3,2);legend('Location','NW');hold on; 
        plot(ax2,h*1e-3,q_rad_TS,'-b','LineWidth',2,'DisplayName','Radiative heat transfer')
        plot(ax2,h*1e-3,q_conv,'-r','LineWidth',2,'DisplayName','Convective heat transfer')
%         ax = gca;
%         ax.YScale = 'log';
        xlabel('Altitude $[km]$')
        ylabel('Heat transfer $[\frac{W}{cm^2}]$')
        grid on;grid minor
        
        ax3 = subplot(2,3,3);legend('Location','NE');hold on; 
        plot(ax3,t,q_rad_TS,'-b','LineWidth',2,'DisplayName','Radiative heat transfer')
        plot(ax3,t,q_conv,'-r','LineWidth',2,'DisplayName','Convective heat transfer')
%         ax = gca;
%         ax.YScale = 'log';
        xlabel('Descent time $[s]$')
        ylabel('Heat transfer $[\frac{W}{cm^2}]$')
        grid on;grid minor
       
        ax4 = subplot(2,3,4);
        plot(ax4,v,T_w,'LineWidth',2)
        xlabel('Velocity $[\frac{km}{s}]$')
        ylabel('Wall temperature $[K]$')
        grid on;grid minor
        
        ax5 = subplot(2,3,5);
        plot(ax5,h*1e-3,T_w,'LineWidth',2)
        xlabel('Altitude $[km]$ ')
        ylabel('Wall temperature $[K]$')
        grid on;grid minor
        
        ax6 = subplot(2,3,6);
        plot(ax6,t,T_w,'LineWidth',2)
        xlabel('Descent time $[s]$ ')
        ylabel('Wall temperature $[K]$')
        grid on;grid minor
        
       
        
end