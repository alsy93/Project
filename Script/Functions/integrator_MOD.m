function [t, y, collectorBank, collectortime] = integrator_MOD(y0,time,par)
%This function solve a reentry dynamic problem by considering a spherical
%coordinate system which pointing the centre of Earth.This is considered as
%an ellipsoid
% INPUTs:
%
%y0                     [1 x 6]                         Initial values of
%                                                       states which are
%                                                       velocity, flight
%                                                       path angle,
%                                                       altitude, latitude,
%                                                       longitude, heading
%                                                       angle ([km/s] [grad]
%                                                       [km] [grad] [grad] 
%                                                       [grad])
%
%time                   [1 x 2]                         Initial and final
%                                                       value of time for
%                                                       simulation. The
%                                                       final time should
%                                                       be when parachute
%                                                       will be opened ([s]
%                                                       [s])
%                                                     
%
%par                    [struct]                        Parameters for the
%                                                       dynamics

% OUTPUTs:
%
%y                      [N x 6]                         State matrix 
%
%t                      [N x 1]                         Time integration
%                                                       vector
%
%collectorBank          [L x 1]                         Controlled variable
%                                                       computed every 4
%                                                       time steps for
%                                                       ode45
%
%collectortime          [L x 1]                         Discrete time
%                                                       vector for bank
%                                                       control

%Check inputs
if nargin ~= 4
    error('Incorrect number of inputs.  See help ecef2enu_12.')
end
if size(v_sim,2) ~= 1
      error('Check the help of this function')
end
if size(hea,2) ~= 1
      error('Check the help of this function')
end
if size(pitch,2) ~= 1
      error('Check the help of this function')
end
if size(bank,2) ~= 1
      error('Check the help of this function')
end

%Checking to see if length of ECEF matrix is the same as the length of the time vector
N = size(v_sim,1);
if N ~= length(hea) || N ~= length(pitch) || N ~= length(bank)
    error('Column size of velocity body vector not equal to size of heading, pitching and banking vectors. Check inputs.')
end


% Inizialize variables for the control of the bank angle

    memOldVel = [];
    initialBank = 60; % [°]
    bank = initialBank;
    memOldBank = initialBank;
    dBank_S0 = -3; % [°] -3
    dBank_S1 = +5; % [°]-4
    eventCount8021_S0 = 0;
    eventCount8021_S1 = 0;
    tempDBank = 0;
    
    collectortime = [];
    collectorBank = [];
 
% System od ODEs
    function dydt = Mechanicalsystm(t,y,par)
            
            v = y(1);
            gamma = y(2);
            h = y(3);
            lat = y(4);
            %long = y(5);
            hea = y(6);
            
            rho = varrho(h);
            
            [Me,Ne] = MeNe(lat,par);

            dydt = zeros(6,1);
            
            dydt(1,1) = -(rho/(2*par.beta))*v^2 - par.g*sind(gamma);
            dydt(2,1) = (((rho/(2*par.alpha))*v*cosd(bank) + cosd(gamma) * (v/(par.Re+h) - par.g/v)))*180/pi;
            dydt(3,1) = v*sind(gamma);
            dydt(4,1) = (v*cosd(gamma)*cosd(hea)/(Me + h))*180/pi;
            dydt(5,1) = (v*cosd(gamma)*sind(hea)/((Ne + h)*cosd(lat)))*180/pi;
            dydt(6,1) = -(((rho/(2*par.alpha)))*v*sind(bank)/cosd(gamma) + ((v/(par.Re+h)))*cosd(gamma)*sind(hea)*tand(lat))*180/pi;
            
    end

    

% Find the corresponding bank angle
        
    function status = myOutputFcn(t,y,flag,par)
        
        switch flag    %No control
            case 'init'
                v = y(1);
                h = y(3);
                
                collectortime = [collectortime;t(1)];
                collectorBank = [collectorBank;initialBank];
            case []    %Control case
                v = y(1);
                h = y(3);
                
                if h >= 80
                    bank = initialBank;
                elseif h>21 && h<80
                    if isempty(memOldVel)
                           memOldVel = v;
                    elseif abs(memOldVel - v) >= 0.2 && abs(memOldVel - v) <= 0.3
                           memOldVel = v;

                        
                            if v <= 4.251 && v >=4.05
                              
                                eventCount8021_S0 = 0;
                                eventCount8021_S1 = eventCount8021_S1 + 1;

                                
                                for i = 1:eventCount8021_S1
                                   tempDBank =  dBank_S1;
                                end

                                S = 1;
                            elseif v < 4.05
                                
                                eventCount8021_S0 = 0;
                                eventCount8021_S1 = eventCount8021_S1 + 1;

                                
                                for i = 1:eventCount8021_S1
                                   tempDBank =  dBank_S1;
                                end
                                
                                tempDBank =  dBank_S1;
                                
                                S = 0;
                            else
                                eventCount8021_S1 = 0;
                                eventCount8021_S0 = eventCount8021_S0 + 1;

                                
                                for i = 1:eventCount8021_S0
                                   tempDBank =  dBank_S0;
                                end

                                S = 0;
                            end

                        
                        bank = (1 - 2*S)*(bank + tempDBank);
                        
                        memOldBank = bank;
                    else
                        bank = memOldBank;
                    end

                elseif h<=21
                    bank = memOldBank;
                end
                
                collectortime = [collectortime;t(end)];
                collectorBank = [collectorBank;memOldBank];
                
            otherwise %Therefore done
        end
        
        status = 0;
        
    end

% Find the corresponding density at each evaluation


    function rho = varrho(h)
           
        if h > 10 
                   h = h*1e3;
              [rho, ~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
              rho = rho*1e9;                     %UM: [kg/km^3]
        else
            h = h*1e3;
            [~, ~, ~, rho] = atmosisa(h);         %UM: [kg/m^3]
            rho = rho*1e9;                      %UM: [kg/km^3]
        end
    end
   

    function [Me,Ne] = MeNe(lat,par)
        
              Me = (par.Re*(1-par.e^2))./((1-par.e^2.*sin(lat).^2).^(3/2));
              Ne = par.Re./sqrt(1-par.e^2.*sin(lat).^2);
    end

% Set options and start modeling simulation


options = odeset('OutputFcn',@myOutputFcn,'RelTol',1e-13,'AbsTol',1e-13,'NormControl','on','Stats','on',...
                 'Events',@h_event);

[t1,y1]= ode45(@Mechanicalsystm,time,y0,options,par);

t = t1;y = y1;


[Me,Ne] = MeNe(y(:,4),par);
rho = varrho (y(:,3));
averBank = mean(collectorBank);


% Plotting of the results

    figure('Name','States as a function of time')
        
        plot(y(:,4),Me./par.Re,y(:,4),Ne./par.Re,'LineWidth',2)
        legend('percentage variation of the meridian radius of curvature ME',...
            'percenage variation of the prime vertical radius of curvature NE')
        xlabel('Latitude $[grad]$');ylabel('Me/Rearth, Ne/Rearth')
        grid on
        grid minor
        hold on
    
    figure('Name','States as a function of time'); 
    
        ax1 = subplot(2,3,1);
        plot(ax1,t,y(:,1))
        title('Velocity')
        grid on;grid minor;hold on
        
        ax2 = subplot(2,3,2);
        plot(ax2,t,y(:,2))
        title('Flight path angle')
        grid on;grid minor;hold on
        
        ax3 = subplot(2,3,3);
        plot(ax3,t,y(:,3))
        title('Altitude')
        grid on;grid minor;hold on
        
        ax4 = subplot(2,3,4);
        plot(ax4,t,y(:,4))
        title('Latitude')
        grid on;grid minor;hold on
        
        ax5 = subplot(2,3,5);
        plot(ax5,t,y(:,5))
        title('Longitude')
        grid on;grid minor;hold on
        
        ax6 = subplot(2,3,6);
        plot(ax6,t,y(:,6))
        title('Heading angle')
        grid on;grid minor;hold on
      
    figure ('Name','Density VS Altitude')  
        plot(rho,y(:,3))
        title('Vehicle deceleration and atmospheric density')
        xlabel('Density $[\frac{kg}{km^3}]$')
        ylabel('Altitude $[km]$')
        grid on;grid minor;hold on
        
    figure ('Name','Density VS time')
        
        plot(t,rho/1e9)
        title('Variation of density in time $[\frac{kg}{m^3}]$')
        grid on;hold on
        
    figure ('Name','Bank angle control')
        
        plot(collectortime,collectorBank,collectortime,averBank)
        title('Variation of bank angle in time')
        grid on; grid minor; hold on
        

end
