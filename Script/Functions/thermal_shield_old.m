function [T_cab,x] = thermal_shield_old(h,v,T_w,q_rad_TS,q_conv,parT,t)

% Change dimension of convective and radiative flux

q_rad_TS = q_rad_TS.*1e4;   %[W/m^2]
q_conv = q_conv.*1e4;       %[W/m^2]

% T_W = T0

%Material response is neglected

%Soyuz TPS characteristics  TO BE CHECKED
% Ablative aeroshell
th_abl = 0.3%;0.028;             % ablation shield thickness: is the X total [m]
qabl_lim = 10000/(1e-4);    % Limit heat flux of aeroshell [W/m^2]
T_m = 800;                  % Melting temperature ablative shield [K]
rho_abl = 1450;             % Density of the ablation material [kg/m^3]
ql_fus = 300*1e3;           % Latent heat of fusion (ablative material) [J/kg]
k_abl = 0.37;               % Thermal conductivity [W/mK]
Cp = 1256.04;               % Specific heat [J/kgK]
alpha = 0.8;                % Solar absorptivity of TPS

% VIM low density insulator
th_ins = 0.008;                         % insulator thickness after TPS
k_ins = 1.5;
Rcond_ins = th_ins/(parT.A_ref*k_ins);

% AlMg6
th_al = 0.002;                        % Alluminimum alloy thickness [m]
k_al = 108;                           % Thermal conductivity of the AlMg [W/mK] 
Rcond_al = th_al/(parT.A_ref*k_al);   % Conductive steel resistance

% Inside the capsule
th_s = 0.03;                                            % Thickness of the structure in steel [m] 
k_s = 50;                                               % Thermal conductivity of steel (W/m-K)
h_in = 20;                                              % Heat transfer coefficient (W/m^2-K)
Rtot_in = 1/(parT.A_ref*h_in) + th_s/(parT.A_ref*k_s);  % Total resistance instance resistance

%=======================================================================================
  % Solve the evolution of temperature inside TPS 
  
    s = abl_rate_consumption(h,v,q_conv,T_w,rho_abl,ql_fus);
    
    function [T_cab,T_endabl, x] = solve_TPS(T_w, q_rad_TS,q_conv,s,parT,t)
    
        L = length(T_w);
        T_cab =  zeros(L,1);
        T_endabl  =  zeros(L,1);
        Rabl =  zeros(L,1);
        Reqtot = zeros(L,1);
        
        x = zeros(L,1);
        x(1) = th_abl;
        
        deltaX = zeros(L,1);
        deltat = zeros(L,1);
        
        Rtot = Rcond_ins + Rcond_al + Rtot_in;
        
            for i = 1:L
                
                if (i > 1)
                    if (T_w(i) >= T_m)
                        
                            deltat(i) = (t(i) - t(i-1));
                            deltaX(i) =  s(i-1)*deltat(i);
                            x(i) = x(i-1) - deltaX(i);
                    else
                        
                            deltaX(i) = 0;
                            x(i) = x(i-1);
                            
                    end
                end
                    
                
                Rabl(i) = x(i)/(parT.A_ref*k_abl);
                Reqtot(i) = (Rabl(i) + Rtot);
                T_cab(i) =  T_w(i) - Reqtot(i)*parT.A_ref*(q_conv(i) + alpha* q_rad_TS(i) - parT.sigma*parT.emiss*T_w(i)^4);% - s(i).*rho_abl.*ql_fus); 
                T_endabl(i) = T_cab(i) + Rtot*(q_conv(i) + alpha* q_rad_TS(i) - parT.sigma*parT.emiss*T_w(i)^4);% - s(i).*rho_abl.*ql_fus);
              
            end
        
    end

   
           
%             for j = 1:L
%                 T_N(j) = T_cab - Rtot*(q_conv(j) + alpha* q_rad_TS(j) - parT.sigma*parT.emiss*T_w(j)^4 + s(j).*rho_abl.*ql_fus);%./kabl.*(th_abl - x(j)); 
%                 q_cab(j) = (T_N(j) - T_cab)/Rtot;
%                
%                 q_cond(j) = q_cab(j);
%               
%                    for i = 1:N
%             
%                        %x(i) = x(i)  - s(j)*(t(j+1)-t(j));  % Position of each node 
%                        T_abl(i,j) = T_N(j) + q_cond(j)./kabl*(th_abl - x(i));
% 
%                         if T_abl(i,j) >= T_m 
% 
%                                 x(i) = x(i) + DELTAx;
%                         else
%                                 x(i) = x(i);
% 
%                         end
%                      thic_plot=[thic_plot;th_abl-x(i)];    
%                    end
%             end
%             
%                 
%        
%          
%       
%          
%         
%         end
       
    
   
% Find corresponding density, speed of sound and temperature of air 
     function [rho, Tatm, a] = varrho(h)

            if h > 10 
                       h = h*1e3;
               	[rho, a, Tatm, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]

            else
                h = h*1e3;
                [Tatm, a, ~, rho] = atmosisa(h);       %UM: [kg/m^3]
            end
     end

% Find corresponding rate of consumption

    function s = abl_rate_consumption(h,v,q_conv,T_w,rho_abl,ql_fus)
        
        [~, Tatm, a] = varrho(h);
        
        Cpatm = 1010;  %Specific heat of air [J/kgK]
        gammaatm = 1.4;
        M = (v.*1e3)./a;
        Htot = Cpatm.*Tatm.*( 1 +(gammaatm-1).*M.^2./2);
        
        q_hw = q_conv.*(1 - Cp.*T_w./Htot);
        
        if sum(q_hw < 0)
           tempIndex = find(q_hw < 0);
           q_hw(tempIndex) = q_conv(tempIndex); 
        end
        
        s = q_hw./(rho_abl*ql_fus);
        
    end

% Plot results
[T_cab,T_endabl, x] = solve_TPS(T_w, q_rad_TS,q_conv,s,parT,t);


figure()
        ax1 = subplot(2,1,1);legend('Location','SE'); hold on; 
        plot(ax1,t,T_cab,'-b','LineWidth',2,'DisplayName','Temperature inside the cabin')
        plot(ax1,t,T_endabl,'-r','LineWidth',2,'DisplayName','Temperature at the end of TPS')
        xlabel('Time $[s]$')
        ylabel('Temperature at the end of ablative coating')
        grid on;grid minor;

        ax2 = subplot(2,1,2);legend('Location','SE'); hold on; 
        plot(ax2,t,x,'LineWidth',2,'DisplayName','Thickness variation in time')
        xlabel('Time $[s]$')
        ylabel('Variation of thicnkess of TPS $[m]$')
        grid on;grid minor;


        %plot(ax1,v,q_conv,'-r','LineWidth',2,'DisplayName','Convectie heat transfer')
%         ax = gca;
%         ax.YScale = 'log';
%         xlabel('Thickness of TPS $[cm]$')
%         ylabel('Temperature $[K]$')
        


end
