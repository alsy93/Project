function [q_conv, T_w]=convective_heat_flux(z,v_AERO,z_i,z_f)
%
% FUNCTION PER TROVARE IL FLUSSO DI CALORE CONVETTIVO NEL PUNTO DI RISTAGNO
% E LA TEMPERATURA DELLA PARETE
%
% INPUT:
% z : altitude
% v : velocità nella terna aerodinamica [3,N]
% z_i : Initial altitude [m]
% z_f : Final altitude [m]
%
% OTPUT:
% q_conv : flusso convettivo nel punto di ristagno
% T_w : wall temperature at the stagnation point
%
% CALLED FUNCTIONS:
% atmosisa
%
% ----------------------------------------

% Costant values
R_n=      % Soyuz nose Radius
K=        % Empirical scaling factor for convective heat flux (see "Simulation of convective heat flux and heat penetration")
h=

% Vectors initialization
n=length(v_AERO);
q_conv=zeros(n,1);
T_w=zeros(n,1);
z_vect=linspace(z_i,z_f,n);

% Classic equation of the convective heat flux
for i=1:n
    %Proprieties of the air, given the altitude
    [T_air(i), ~, ~, rho(i)] = atmosisa(z_vect(i))   %T_air[K] -  rho[g/m^3]
    
    % Heat flux at the stagnation point at that altitude
    q_conv(i)=K*(rho(i)/R_n)^(1/2)*v(i);
    
    % Wall temperature at the stagnation point at that altitude
    T_w(i)=T_air(i)+q_conv(i)/h;
end

%%% PROBLEMI
% 1 - Devo trovare come scrivere h (convective heat flux), perché il valore
%     non è costante e dipende da velocità, temperatura, tipo di fluido,
%     cazzi e ammazzi e quindi devo trovare una correlazione adatta.
% 2 - Poi mi sono accorta che la function atmosisa è limitata alla
%     troposfera e la nostra analisi comincia ben più in alto quindi devo
%     trovare uno (o più) altri modelli da cui recuperare temperatura e
%     densità dell'aria
    
    
    
    
    
    
    
    
    
    
    
    