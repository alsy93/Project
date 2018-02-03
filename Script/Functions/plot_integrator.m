function [vargout] = plot_integrator(t45,y45,bank45,timebank45,t113,y113,bank113,timebank113)
    

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

figure ('Name','Step size variation in integration time with a set Max step size 0.3')
        
        legend();hold on
        plot(t45,deltaStep45,'LineWidth',2,'DisplayName','Default step size variation of ODE45 ')
        plot(t113,deltaStep113,'LineWidth',2,'DisplayName','Default step size variation of ODE113 ')
        grid on; grid minor
        
figure('Name','States as a function of time of ODE45');
       
        ax1 = subplot(2,3,1);
        plot(ax1,t45,y45(:,1),'LineWidth',2,'DisplayName','ODE45');
         
        title('Velocity')
        grid on;grid minor;hold on
        
        ax2 = subplot(2,3,2);
        plot(ax2,t45,y45(:,2),'LineWidth',2);
        title('Flight path angle')
        grid on;grid minor;hold on
        
        ax3 = subplot(2,3,3);
        plot(ax3,t45,y45(:,3),'LineWidth',2);
        title('Altitude')
        grid on;grid minor;hold on
        
        ax4 = subplot(2,3,4);
        plot(ax4,t45,y45(:,4),'LineWidth',2);
        title('Latitude')
        grid on;grid minor;hold on
        
        ax5 = subplot(2,3,5);
        plot(ax5,t45,y45(:,5),'LineWidth',2);
        title('Longitude')
        grid on;grid minor;hold on
        
        ax6 = subplot(2,3,6);
        plot(ax6,t45,y45(:,6),'LineWidth',2);
        title('Heading angle')
        grid on;grid minor;hold on
        
        
   vargout(2) = figure('Name','States as a function of time of ODE113');
       
        ax1 = subplot(2,3,1);
        plot(ax1,t113,y113(:,1),'LineWidth',2,'DisplayName','ODE45');
         
        title('Velocity')
        grid on;grid minor;hold on
        
        ax2 = subplot(2,3,2);
        plot(ax2,t113,y113(:,2),'LineWidth',2);
        title('Flight path angle')
        grid on;grid minor;hold on
        
        ax3 = subplot(2,3,3);
        plot(ax3,t113,y113(:,3),'LineWidth',2);
        title('Altitude')
        grid on;grid minor;hold on
        
        ax4 = subplot(2,3,4);
        plot(ax4,t113,y113(:,4),'LineWidth',2);
        title('Latitude')
        grid on;grid minor;hold on
        
        ax5 = subplot(2,3,5);
        plot(ax5,t113,y113(:,5),'LineWidth',2);
        title('Longitude')
        grid on;grid minor;hold on
        
        ax6 = subplot(2,3,6);
        plot(ax6,t113,y113(:,6),'LineWidth',2);
        title('Heading angle')
        grid on;grid minor;hold on
        
            
     
   vargout(3) = figure ('Name','Bank angle control ODE45');
        
        legend(); hold on
        plot(timebank45,bank45,'LineWidth',2)
        title('Variation of bank angle in time of ODE45');
        grid on; grid minor; hold on
 
    vargout(4) = figure ('Name','Bank angle control ODE113');
        
        legend(); hold on        
        plot(timebank113,bank113,'LineWidth',2);
        title('Variation of bank angle in time of ODE113');
        grid on; grid minor; hold on
end