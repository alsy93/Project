function [vargout] = plot_integrator_0(t45,y45,bank45,timebank45,t113,y113,bank113,timebank113)
    

% For ODE 45

N45 = length(t45);
deltaStep45 = zeros(N45,1);
deltaStep45(1) = t45(1);
for j = 2:N45
    
    deltaStep45(j) =(t45(j) - t45(j-1));
    
end
% For ODE 113

N113 = length(t113);
deltaStep113 = zeros(N113,1);
deltaStep113(1) = t113(1);
for j = 2:N113
    
    deltaStep113(j) =(t113(j) - t113(j-1));
    
end

%Plotting

figure ('Name','Step size variation')
        
        legend();hold on
        vargout(1) = plot(t45,deltaStep45,'LineWidth',2,'DisplayName','Default step size variation of ODE45 ');
        vargout(2) = plot(t113,deltaStep113,'LineWidth',2,'DisplayName','Default step size variation of ODE113 ');
        title('Step size variation in integration time with a default Max step size');
        xlabel('Integration time $[s]$'); ylabel('Value of time step variation $[s]$');
        grid on; grid minor


figure('Name','States of ODE113');
       
        ax1 = subplot(2,3,1);
        vargout(3) = plot(ax1,t113,y113(:,1),'LineWidth',2);
        title('Velocity $[\frac{kg}{km^3}]$')
        xlabel('Time $[s]$');
        grid on;grid minor;hold on
        
        ax2 = subplot(2,3,2);
        vargout(4) = plot(ax2,t113,y113(:,2),'LineWidth',2);
        title('Flight path angle $[grad]$')
        xlabel('Time $[s]$');
        grid on;grid minor;hold on
        
        ax3 = subplot(2,3,3);
        vargout(5) = plot(ax3,t113,y113(:,3),'LineWidth',2);
        title('Altitude$[km]$')
        xlabel('Time $[s]$');
        grid on;grid minor;hold on
        
        ax4 = subplot(2,3,4);
        vargout(6) = plot(ax4,t113,y113(:,4),'LineWidth',2);
        title('Latitude $[grad]$')
        xlabel('Time $[s]$');
        grid on;grid minor;hold on
        
        ax5 = subplot(2,3,5);
        vargout(7) = plot(ax5,t113,y113(:,5),'LineWidth',2);
        title('Longitude $[grad]$')
        xlabel('Time $[s]$');
        grid on;grid minor;hold on
        
        ax6 = subplot(2,3,6);
        vargout(8) = plot(ax6,t113,y113(:,6),'LineWidth',2);
        title('Heading angle $[grad]$')
        xlabel('Time $[s]$');
        grid on;grid minor;hold on
        
            
 figure ('Name','Bank angle ODE45 - ODE113');
        
        legend('Position','NE'); hold on;
        vargout(9) = plot(timebank45,bank45,'LineWidth',2,'DisplayName','Bank angle ODE45');
        vargout(10) = plot(timebank113,bank113,'LineWidth',2,'DisplayName','Bank angle ODE113');
        xlabel('Integration time $[s]$'); ylabel('Bank angle $[grad]$');
        grid on; grid minor; hold on
      
end