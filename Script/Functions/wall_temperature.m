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
% T_w = zeros(1,N);
% q_rad_TS = zeros(1,N);
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
        if h > 11 %after tropopause
              h = h*1e3;
              [rho,~, ~, ~, ~, ~] = atmos(h);   %UM: [kg/m^3]
              rho = rho*1e9;                    %UM: [kg/km^3]
        else 
              h = h*1e3;
              [rho,~, ~, ~,~]=tropos(h);        %UM: [kg/m^3]
              rho = rho*1e9;                    %UM: [kg/km^3]
        end
    end

end