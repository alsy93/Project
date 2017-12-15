function [T_w, q_rad_TS] = wall_temperature(v,h,parT)

% Function to calculate the wall temperature at the stagnation point.
% It comes from an equality of two equations, the analytic Stefan-Boltzman
%equation for radiative heat flux and the statistical Sutton-Tauber one.


tic
rho = varrho(h);

q_rad_TS = (parT.Rn^parT.a) .* (rho.^parT.b) .*(v.^(parT.d));
T_w = (q_rad_TS ./ (parT.sigma * parT.emiss)).^(1/4);


% % Vectors initialization
% N = length(v);
% T_w = zeros(N,1);
% q_rad_TS = zeros(N,1);
% rho = varrho(h);
% 
% %
% for i=1:N
%     q_rad_TS(i) = parT.C * parT.Rn^parT.a * rho(i)^parT.b * v(i)^parT.d;
%     T_w(i) = (q_rad_TS(i) / (parT.sigma * parT.emiss))^(1/4);
% end
toc

% Find the corresponding density at each evaluation

    function rho = varrho(h)
        if h > 11e3 %after tropopause
              [rho,~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
        else 
              [rho,~, ~, ~,~]=tropos(h);        %UM: [kg/m^3]
       
        end
    end
% Plotting  of results

figure(7)
        ax1 = subplot(2,2,1);
        plot(ax1,v,q_rad_TS)
        xlabel('Velocity [km/s]')
        ylabel('Radiative heat transfer [w/cm^2]')
        grid on;grid minor
        
        ax2 = subplot(2,2,2);
        plot(ax2,h*1e-3,T_w)
        xlabel('Altitude [km] ')
        ylabel('Wall temperature [k]')
        grid on;grid minor
        
        ax3 = subplot(2,2,3);
        plot(ax3,v,T_w)
        xlabel('Velocity [km/s]')
        ylabel('Wall temperature [k]')
        grid on;grid minor
        
        ax4 = subplot(2,2,4);
        plot(ax4,h*1e-3,q_rad_TS)
        xlabel('Altitude [km]')
        ylabel('Radiative heat transfer [w/cm^2]')
        grid on;grid minor
        
end