function [vargout] = flight_envelope(h,v,T_w,par,parT)

%Exponentialatmospheric model: rho = rho0*exp(-beta*h) ;
beta = 835/par.Re;      % [1/km]
rho0 = 1.225*1e9;       % [kg/km^3]
% This function provides the flight envelope constraints

% 1)Equilibrium trajectory [gamma is zero and also the dgamma/dt is zero]:
%   maximum lift ceiling. It should be upper the computed solution

    h_eq = log((rho0.*v.^2)./(2*par.alpha.*(par.g - v.^2./par.Re)))./...
           (beta - 1/par.Re .* (2*par.g - (v.^2./par.Re)./(par.g - ...
           (v.^2./par.Re))));
       
% 2)Thermal Limit (In Flux): a maximal
%   value of the reference flux is usually used, i.e. the flux, which a 
%   sphere whose radius is equal to the radius of curvature that the “nose”
%   of the vehicle (Rn).It should be below the computed solution if
%   ablative TPS is not used

    
    T_max = max(T_w);
    
    h_tl = -2/beta.*log((parT.sigma*1e-6*parT.emiss.*T_max.^4*...
           sqrt(parT.Rn*1e-5))./(parT.Ks*sqrt(rho0).*v.^3));

    % By using information about radiative flux, values with v<= 2 km/s are neglected 
    n = length(h_tl);
    for i = 1:n
        if h_tl(i) <= 0
           h_tl(i) = NaN;
        else
            h_tl(i) = h_tl(i);
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

rho = varrho(h);

    n_tot = 0.5.*rho.*par.S*sqrt(par.cl^2+par.cd^2).*v.^2./(par.m*par.g);
    h_lfl = -log(n_tot.*par.g.*par.m./(0.5.*rho0.*v.^2.*par.S.*par.cl))./beta;


% Plot solutions




    figure('Name','Flight envelope');hold on; grid on; grid minor
    title('System and operational constraints');
    xlabel('Relative velocity $[\frac{km}{s}]$');
    ylabel('Altitude $[km]$');
    legend();
    vargout(1) = plot (v,h,'-b','LineWidth',2,'DisplayName','Computed solution');
    vargout(2) = plot(v,h_eq,'-r','LineWidth',2,'DisplayName','Maximum lift ceiling');
               
    %The requirement of ablative TPS is highlighted
    
    N = length(h);
    
    for i = 1:N
        if h_tl(i) >= h(i)
            hTPS = h_tl(i);     %Mark when TPS is required
            vTPS = v(i); 
        end
    end
    
    indexH = find(h_tl<hTPS,1);
    indexV = find(v<vTPS,1);
    
    vargout(3) = plot(v(1:indexV),h_tl(1:indexH),'--g','LineWidth',2,'DisplayName','Thermal limit flux: ablative TPS is required');
    vargout(4) = plot(v(indexV:end),h_tl(indexH:end),'-g','LineWidth',2,'DisplayName','Thermal limit flux: reusable TPS is preferable');
    vargout(5) = plot(vTPS,hTPS,'or','LineWidth',2,'DisplayName','TPS is needed for higher values');
    vargout(6) = plot(v,h_lfl,'-m','LineWidth',2,'DisplayName','Load factor limit');
end



    
    

