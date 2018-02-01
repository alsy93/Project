function q_conv = convective_flux(v,h,parT)
%Fucntion to calculate the convective heat transfer at the stagnation point
%by using the Sutton-Graves theory, which is coming from Chapman equation

% Find corresponding density 
    function rho = varrho(h)
           
        if h > 10000
           [rho, ~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
           
        else
           [~, ~, ~, rho] = atmosisa(h);      %UM: [kg/m^3]
           
        end
    end
rho = varrho(h);
q_conv = (parT.Ks.*(rho./(parT.Rn*1e-2)).^(1/2).*((v.*1e3).^3)).*1e-4;
            
end
