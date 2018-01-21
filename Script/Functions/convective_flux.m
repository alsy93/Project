function q_cond = convective_flux(v,h,parT)
%Fucntion to calculate the convective heat transfer at the stagnation point
%by using the Sutton-Graves theory, which is coming from Chapman equation
tic
% Find corresponding density 
    function rho = varrho(h)
           
        if h > 10000
           [rho, ~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
           
        else
           [~, ~, ~, rho] = atmosisa(h);      %UM: [kg/m^3]
           
        end
    end
rho = varrho(h);
q_cond = parT.Ks.*(rho./parT.Rn).^(1/2).*(v.^3);
toc

% Plotting  of results

figure(9)
        ax1 = subplot(1,2,1);
        plot(ax1,v,q_cond)
        xlabel('Velocity (km/s)')
        ylabel('Convectie heat transfer (w/cm^2)')
        grid on;grid minor
        
        ax2 = subplot(1,2,2);
        plot(ax2,h*1e-3,q_cond)
        xlabel('Altitude (km) ')
        ylabel('Convective heat transfer (w/cm^2)')
        grid on;grid minor
        
       
end
