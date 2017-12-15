function [T_w, q_rad_TS] = wall_temperature(v, h, parT)

% Function to calculate the wall temperature at the stagnation point.
%
% !!! NON ORTHODOX WAY
% We correlate two equations, the analytic Stefan-Boltzman equation for
% radiative heat flux and the statistical Sutton-Tauber one.
%

tic

rho = varrho(h);

q_rad_TS = parT.C .* parT.Rn^parT.a .* rho.^parT.b .* v.^parT.d;
T_w = (q_rad_TS ./ (parT.sigma * parT.emiss)).^(1/4);

% % Vectors initialization
% N = length(v);
% T_w = zeros(N,1);
% q_rad_TS = zeros(N,1);
% rho = varrho(h);
% 
% %
% for i=1:N
%     q_rad_TS(i) = parT.C * parT.Rn^parT.a * rho(i)^parT.b * v(i)^parT.d;
%     T_w(i) = (q_rad_TS(i) / (parT.sigma * parT.emiss))^(1/4);
% end
toc

% Find the corresponding density at each evaluation

    function rho = varrho(h)
        if h > 11000 %after tropopause
              [rho,~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
        else 
              [rho,~, ~, ~,~]=tropos(h);        %UM: [kg/m^3]
        end
    end

end