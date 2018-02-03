function [vargout] = flight_envelope_sensitivity(h1,v1,T_w1,h2,v2,T_w2,par,parT)
%TO BE USED ONLY FOR SENSITIVITY ANALYSIS

%Exponentialatmospheric model: rho = rho0*exp(-beta*h) ;
beta = 835/par.Re;      % [1/km]
rho0 = 1.225*1e9;       % [kg/km^3]
% This function provides the flight envelope constraints

% 1)Equilibrium trajectory [gamma is zero and also the dgamma/dt is zero]:
%   maximum lift ceiling. It should be upper the computed solution

    h_eq_1 = log((rho0.*v1.^2)./(2*par.alpha.*(par.g - v1.^2./par.Re)))./...
           (beta - 1/par.Re .* (2*par.g - (v1.^2./par.Re)./(par.g - ...
           (v1.^2./par.Re))));
       
     
    h_eq_2 = log((rho0.*v2.^2)./(2*par.alpha.*(par.g - v2.^2./par.Re)))./...
           (beta - 1/par.Re .* (2*par.g - (v2.^2./par.Re)./(par.g - ...
           (v2.^2./par.Re))));  
       
% 2)Thermal Limit (In Flux): a maximal
%   value of the reference flux is usually used, i.e. the flux, which a 
%   sphere whose radius is equal to the radius of curvature that the “nose”
%   of the vehicle (Rn).It should be below the computed solution if
%   ablative TPS is not used

    
    T_max_1 = max(T_w1);
    
    h_tl_1 = -2/beta.*log((parT.sigma*1e-6*parT.emiss.*T_max_1.^4*...
           sqrt(parT.Rn*1e-5))./(parT.Ks*sqrt(rho0).*v1.^3));
       
   T_max_2 = max(T_w2);
    
    h_tl_2 = -2/beta.*log((parT.sigma*1e-6*parT.emiss.*T_max_2.^4*...
           sqrt(parT.Rn*1e-5))./(parT.Ks*sqrt(rho0).*v2.^3));


    % By using information about radiative flux, values with v<= 2 km/s are neglected 
    n1 = length(h_tl_1);
    for i = 1:n1
        if h_tl_1(i) <= 0
           h_tl_1(i) = NaN;
        else
            h_tl_1(i) = h_tl_1(i);
        end
    end
    
     % By using information about radiative flux, values with v<= 2 km/s are neglected 
    n2 = length(h_tl_2);
    for i = 1:n2
        if h_tl_2(i) <= 0
           h_tl_2(i) = NaN;
        else
            h_tl_2(i) = h_tl_2(i);
        end
    end
    
% 3)Load factor limit: the maximal load factor depends on the size and the
%   presence of the crew
    %Find corresponding density for a given altitude
     function rho = varrho(h)

            if h > 10 
                       h = h*1e3;
                  [rho, ~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
                  rho = rho*1e9;                     %UM: [kg/km^3]
            else
                h = h*1e3;
                [~, rho] = atmosisa(h);             %UM: [kg/m^3]
                rho = rho*1e9;                      %UM: [kg/km^3]
            end
     end

rho1 = varrho(h1);

    n_tot_1 = 0.5.*rho1.*par.S*sqrt(par.cl^2+par.cd^2).*v1.^2./(par.m*par.g);
    h_lfl_1 = -log(n_tot_1.*par.g.*par.m./(0.5.*rho0.*v1.^2.*par.S.*par.cl))./beta;

    
rho2 = varrho(h2);

    n_tot_2 = 0.5.*rho2.*par.S*sqrt(par.cl^2+par.cd^2).*v2.^2./(par.m*par.g);
    h_lfl_2 = -log(n_tot_2.*par.g.*par.m./(0.5.*rho0.*v2.^2.*par.S.*par.cl))./beta;


% Plot solutions




    figure('Name','Flight envelope');
    
    ax1 = subplot(1,2,1);legend('Location','SE');hold on;
   
    vargout(1) = plot (ax1,v1,h1,'-b','LineWidth',2,'DisplayName','Computed solution');
    vargout(2) = plot(ax1,v1,h_eq_1,'-r','LineWidth',2,'DisplayName','Maximum lift ceiling');
               
    %The requirement of ablative TPS is highlighted
    
    N1 = length(h1);
    
    for i = 1:N1
        if h_tl_1(i) >= h1(i)
            hTPS1 = h_tl_1(i);     %Mark when TPS is required
            vTPS1 = v1(i); 
        end
    end
    
    indexH1 = find(h_tl_1<hTPS1,1);
    indexV1 = find(v1<vTPS1,1);
    
    vargout(3) = plot(ax1,v1(1:indexV1),h_tl_1(1:indexH1),'--g','LineWidth',2,'DisplayName','Thermal limit flux: ablative TPS is required');
    vargout(4) = plot(ax1,v1(indexV1:end),h_tl_1(indexH1:end),'-g','LineWidth',2,'DisplayName','Thermal limit flux: reusable TPS is preferable');
    vargout(5) = plot(ax1,vTPS1,hTPS1,'or','LineWidth',2,'DisplayName','TPS is needed for higher values');
    vargout(6) = plot(ax1,v1,h_lfl_1,'-m','LineWidth',2,'DisplayName','Load factor limit');
    title('System and operational constraints with $\gamma = -0.75$');
    xlabel('Relative velocity $[\frac{km}{s}]$');
    ylabel('Altitude $[km]$');
    grid on; grid minor

    ax2 = subplot(1,2,2);legend('Location','SE');hold on;
   
    vargout(1) = plot(ax2,v2,h2,'-b','LineWidth',2,'DisplayName','Computed solution');
    vargout(2) = plot(ax2,v2,h_eq_2,'-r','LineWidth',2,'DisplayName','Maximum lift ceiling');

        %The requirement of ablative TPS is highlighted

        N2 = length(h2);

        for i = 1:N2
            if h_tl_2(i) >= h2(i)
                hTPS2 = h_tl_2(i);     %Mark when TPS is required
                vTPS2 = v2(i); 
            end
        end

        indexH2 = find(h_tl_2<hTPS2,1);
        indexV2 = find(v2<vTPS2,1);

         
    vargout(3) = plot(ax2,v2(1:indexV2),h_tl_2(1:indexH2),'--g','LineWidth',2,'DisplayName','Thermal limit flux: ablative TPS is required');
    vargout(4) =plot(ax2,v2(indexV2:end),h_tl_2(indexH2:end),'-g','LineWidth',2,'DisplayName','Thermal limit flux: reusable TPS is preferable');
    vargout(5) =plot(ax2,vTPS2,hTPS2,'or','LineWidth',2,'DisplayName','TPS is needed for higher values');
    vargout(6) =plot(ax2,v2,h_lfl_2,'-m','LineWidth',2,'DisplayName','Load factor limit');
    title('System and operational constraints $\gamma = -3.5$');
    xlabel('Relative velocity $[\frac{km}{s}]$');
    ylabel('Altitude $[km]$');
    grid on; grid minor

end

