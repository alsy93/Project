function [t113_1, y113_1,bank113,timebank113] = integrator_MOD113_1(y0,time,par)
tic
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


%2) With a set Max step size
options1 = odeset('OutputFcn',@myOutputFcn,'RelTol',1e-13,'AbsTol',1e-13,'NormControl','on','Stats','on',...
                 'Events',@h_event,'MaxStep',0.3);             
         
disp('=================================================================='); 
disp('ODE113 evaluations with a Max step set at 0.3');             
[t113_1,y113_1]= ode113(@Mechanicalsystm,time,y0,options1,par);
bank113 = collectorBank;
timebank113 = collectortime;
Integrator_timeODE113 = toc;
fprintf('Run time of integration with a Max step size 0.3 is %d.\n',Integrator_timeODE113);



end