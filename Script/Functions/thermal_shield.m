function [T_cab,x] = thermal_shield(h,v,T_w,q_rad_TS,q_conv,parT,t)

%Note that q_rad_Ts and q_conv should be in W/m^2

%Soyuz TPS characteristics 

% Ablative aeroshell
th_abl = 0.028;             % ablation shield thickness: is the X total [m]
T_m = 755;                  % Melting temperature ablative shield [K]
rho_abl = 1450;             % Density of the ablation material [kg/m^3]
ql_fus =29600*1e3;          % Latent heat of fusion (ablative material) [J/kg]
k_abl = 0.37;               % Thermal conductivity [W/mK]
Cp = 1256.04;               % Specific heat [J/kgK]
alpha = 0.9;                % Solar absorptivity 

% VIM low density insulator
th_ins = 0.008;                         % insulator thickness after TPS
k_ins = 0.15;                           % Thermal conductivity [W/mK]
Rcond_ins = th_ins/(parT.A_ref*k_ins);  % Conductive resistance of insulator
rho_ins = 207;                          % Density [kg/m^3]
Cp_ins = 1590;                          % Specific heat [J/kgK]
lm_ins = rho_ins*parT.A_ref*th_ins*Cp_ins;

% AlMg6
th_al = 0.002;                        % Alluminimum alloy thickness [m]
k_al = 108;                           % Thermal conductivity of the AlMg [W/mK] 
Rcond_al = th_al/(parT.A_ref*k_al);   % Conductive Al-Mg resistance 
Cp_al = 1040;                         % Specific heat [J/kgK]
rho_al = 1900;                        % Density [kg/m^3]
lm_al = rho_al*parT.A_ref*th_al*Cp_al;


% Soyuz structural and internal characteristics
th_steel = 0.03;                            % Thickness of the structure in steel [m] 
k_s = 50;                                   % Thermal conductivity of steel (W/mK)
Rcond_steel=  th_steel/(parT.A_ref*k_s);    % Conductive steel resistance 
Cp_steel = 490;                             % Specific heat [J/kgK]
rho_steel = 8050;                           % Density [kg/m^3]
lm_steel = rho_steel*parT.A_ref*th_steel*Cp_steel;

h_cab = 20;                           % Heat transfer coefficient (W/m^2K)
Rconv_cab = 1/(parT.A_ref*h_cab);      % Convective air resistance inside cabin 
Cp_air= 1010;                         % Specific heat of air [J/kgK]
rho_air = 1.225;                      % Density [kg/m^3]
lm_cab = rho_air*parT.Vol*Cp_air;    

%=======================================================================================
  % Solve the evolution of temperature inside TPS 

  s = abl_rate_consumption(h,v,q_conv,T_w,rho_abl,ql_fus);
    
 % Set nodes temperatures (Guessed values)
 T_N2 = 700;        %Lumped inside ablative material [K]
 T_N3 = 600;        %Lumped contact end of ablative material and insulator [K]
 T_N4 = 400;        %Lumped contact insulator and Al-Mg alloy substrate[K]
 T_N5 = 350;        %Lumped contact Al substrate-steel structure[K]
 T_N6 = 349;        %Lumped structure internal wall [K]
 T_N7 = 345; 
 T_N8 = 340;
 T_N9 = 330;
 T_cab = 298;       %Initial cabin temperature [K]
  
 L = length(t); 
 T_nodes0 = [T_w(1) T_N2 T_N3 T_N4 T_N5 T_N6 T_N7 T_N8 T_N9 T_cab];
 T_nodes = zeros(L,10);
 T_nodes(1,:) = T_nodes0;
 
 Qin = zeros(L,1);
 Qin(1) = parT.A_ref*(q_conv(1) +  q_rad_TS(1)); 
 
 deltat = zeros(L,1);
 deltaX = zeros(L,1);
 x = zeros(L,1);
 x(1) = th_abl;
 
 R_abl = zeros(L,1);
 Rc_ablis = zeros(L,1);

% Contact resistance

Rc_insAl = (Rcond_ins + Rcond_al)/4;
Rc_AlSt = (Rcond_al + Rcond_steel)/4;
Rc_StCab = (Rcond_steel + Rconv_cab)/4;
 
 %Compute values for each i-time discretization
 
 for i = 2:L
     
      deltat(i) = t(i) - t(i-1);
      Qin(i)= parT.A_ref*(q_conv(i) +  q_rad_TS(i));% - parT.sigma*parT.emiss*T_w(i)^4 - s(i).*rho_abl.*ql_fus);
      
 %ABLATIVE     
 %Thickness variation for pyrolysis
  
    if T_w(i) >= T_m
        
        deltaX(i) =  s(i-1)*deltat(i);
        x(i) = x(i-1) - deltaX(i);
    else
        deltaX(i) = 0;
        x(i) = x(i-1);
 
    end
 %Find corresponding ablative equivalence resistance  
 R_abl(i) = x(i)/(2*(parT.A_ref*k_abl));
    
 %Lumped mass 
 lm_abl = (rho_abl*parT.A_ref*x(i))*Cp;
 
 % In the node #1
 T_nodes(i,1) = T_w(i);
 
 % In the node #2
 Rc_ablis(i) = (R_abl(i) + Rcond_ins/2) /2;
 T_nodes(i,2) = T_nodes(i-1,2) + deltat(i)/(lm_abl) * ((T_w(i-1) - T_nodes(i-1,2))/R_abl(i) - (T_nodes(i-1,2) - T_nodes(i-1,3))/(R_abl(i) + Rc_ablis(i)) - s(i-1)*rho_abl*ql_fus );
 
 % In the node #3

 T_nodes(i,3) = T_nodes(i-1,3) + 2*deltat(i)/(lm_abl +lm_ins) * ((T_nodes(i-1,2) - T_nodes(i-1,3))/(R_abl(i) + Rc_ablis(i)) - (T_nodes(i-1,3) - T_nodes(i-1,4))/(Rc_ablis(i) + Rcond_ins/2) );
 
 %INSULATOR
 
 % In the node #4
 T_nodes(i,4) = T_nodes(i-1,4) + deltat(i)/(lm_ins) * ((T_nodes(i-1,3) - T_nodes(i-1,4))/(Rc_ablis(i) + Rcond_ins/2) - (T_nodes(i-1,4) - T_nodes(i-1,5))/(Rcond_ins/2 + Rc_insAl) );
 
 % In the node #5
 T_nodes(i,5) = T_nodes(i-1,5) + 2*deltat(i)/(lm_ins + lm_al) * ((T_nodes(i-1,4) - T_nodes(i-1,5))/(Rcond_ins/2 + Rc_insAl) - (T_nodes(i-1,5) - T_nodes(i-1,6))/(Rcond_ins/2 + Rc_insAl) );
 
 %Al_mg ALLOY
 
 T_nodes(i,6) = T_nodes(i-1,6) + deltat(i)/(lm_al) * ((T_nodes(i-1,5) - T_nodes(i-1,6))/(Rcond_ins/2 + Rc_insAl) - (T_nodes(i-1,6) - T_nodes(i-1,7))/(Rcond_al/2 + Rc_AlSt));

 
 T_nodes(i,7) = T_nodes(i-1,7) + 2*deltat(i)/(lm_al + lm_steel) * ((T_nodes(i-1,6) - T_nodes(i-1,7))/(Rcond_al/2 + Rc_AlSt) - (T_nodes(i-1,7) - T_nodes(i-1,8))/(Rc_AlSt + Rcond_steel/2) );
    
 % STEEL 
 T_nodes(i,8) = T_nodes(i-1,8) + deltat(i)/(lm_steel) * ((T_nodes(i-1,7) - T_nodes(i-1,8))/(Rc_AlSt + Rcond_steel/2) - (T_nodes(i-1,8) - T_nodes(i-1,9))/(Rcond_steel/2 + Rc_StCab) );  
 
 T_nodes(i,9) = T_nodes(i-1,9) + 2*deltat(i)/(lm_steel + lm_cab) * ((T_nodes(i-1,8) - T_nodes(i-1,9))/(Rcond_steel/2 + Rc_StCab) - (T_nodes(i-1,9) - T_nodes(i-1,10))/(Rconv_cab + Rc_StCab));
 
 %CABIN
  
 T_nodes(i,10) = T_nodes(i-1,10) + deltat(i)/(lm_cab) *((T_nodes(i-1,9) - T_nodes(i-1,10))/(Rconv_cab + Rc_StCab) );
 end
  
  
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


%===========================================================================

T_abl1 = T_nodes(:,2);
T_abl2 = T_nodes(:,3);
T_insAl = T_nodes(:,6);
T_AlSt = T_nodes(:,8);
T_cab = T_nodes(:,10);

[T_cabmax, indexTcabinmax] = max(T_cab);
altitude_Tmax = h(indexTcabinmax);
fprintf('Altitude for which the temperature inside the cabin is maximum is %f km. \n',altitude_Tmax)
fprintf('The maximum temperature inside the cabin is %f K.\n',T_cabmax);

%Plot results

figure('Name','TPS temperatures')
        ax1 = subplot(2,1,1);legend('Location','NE'); hold on; 
        plot(ax1,t,T_cab,'-b','LineWidth',2,'DisplayName','Temperature inside the cabin')
        plot(ax1,t,T_abl1,'-r','LineWidth',2,'DisplayName','Temperature inside TPS at node 2')
        plot(ax1,t,T_abl2,'--r','LineWidth',2,'DisplayName','Temperature at the end of TPS')
        xlabel('Time $[s]$')
        ylabel('Temperature $[K]$')
        grid on;grid minor;

        ax2 = subplot(2,1,2);legend('Location','NE'); hold on; 
        plot(ax2,t,x,'LineWidth',2,'DisplayName','Thickness variation in time')
        xlabel('Time $[s]$')
        ylabel('Variation of thicnkess of TPS $[m]$')
        grid on;grid minor;
        
figure('Name','Temperature in thickness')
        
        totalTK =[0, x(end)/2, x(end), x(end)+th_ins/2, x(end)+th_ins, x(end)+th_ins+th_al/2, x(end)+th_ins+th_al, x(end)+th_ins+th_al+th_steel/2, x(end)+th_ins+th_al+th_steel, x(end)+th_ins+th_al+th_steel ]; 
        
        legend('Location','NE'); hold on;
        plot(totalTK,T_nodes(indexTcabinmax,:),'LineWidth',2,'DisplayName','Temperature variation inside whole thickness');
        title('Temperature at nodes when in the cabin is max');
        xlabel('Thickness $[m]$')
        ylabel('Temperature $[K]$')
        grid on; grid minor


end