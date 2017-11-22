function [q_conv, T_w]=convective_heat_flux(z,v_AERO,z_i,z_f)
%
% FUNCTION PER TROVARE IL FLUSSO DI CALORE CONVETTIVO NEL PUNTO DI RISTAGNO
% E LA TEMPERATURA DELLA PARETE
%
% INPUT:
% z : altitude
% v : velocità nella terna aerodinamica [3 x N]
% z_i : Initial altitude [m]
% z_f : Final altitude [m]
%
% OTPUT:
% q_conv : flusso convettivo nel punto di ristagno [1 x N]
% T_w : wall temperature at the stagnation point   [1 x N]
%
% CALLED FUNCTIONS:
% atmosisa
%atmosnrlmsise00
%
% ----------------------------------------

% Costant values
R_n=      % Soyuz nose Radius
K=        % Empirical scaling factor for convective heat flux (see "Simulation of convective heat flux and heat penetration")
h=

% Vectors initialization
N=length(v_AERO);
q_conv=zeros(1,N);
T_w=zeros(1,N);
T_air=zeros(1,N);
rho=zeros(1,N);
%z_vect=linspace(z_i,z_f,N);

% Classic equation of the convective heat flux
for i=1:N
    %Proprieties of the air, given the altitude
    [T_air(1,i), ~, ~,rho(1,i)] = atmosisa(norm(z(:,i)))   %T_air[K] -  rho[g/m^3]
    [T, rho_] = atmosnrlmsise00(h,lat,lon,year,doy,sec,varargin)
    % Heat flux at the stagnation point at that altitude
    q_conv(1,i)=K*(rho(1,i)/R_n)^(1/2)*norm(v_AERO(:,i));
    
    % Wall temperature at the stagnation point at that altitude
    T_w(1,i)=T_air(1,i)+q_conv(1,i)/h;
end

%%% PROBLEMI
% 1 - Devo trovare come scrivere h (convective heat flux), perché il valore
%     non è costante e dipende da velocità, temperatura, tipo di fluido,
%     cazzi e ammazzi e quindi devo trovare una correlazione adatta.
% 2 - Poi mi sono accorta che la function atmosisa è limitata alla
%     troposfera e la nostra analisi comincia ben più in alto quindi devo
%     trovare uno (o più) altri modelli da cui recuperare temperatura e
%     densità dell'aria
% bisogna trattare bene quota perchè basta z come vettore senza valore
% iniziale e finale che sono già presenti nel medesimo
    
    
    
    
    
    
    
    
    
    
    