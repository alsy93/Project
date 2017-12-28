function rho = varrho(h,lat,long)
       
        if h >= 86
             h = h*1e3;
             t_0 = [2017 06 02 16 47 25]; 
             timezone = +3;
             UTC = t_0(6) + t_0(5)*60 + (timezone + t_0(4))*3600 + ... 
                   t_0(3)*24*3600 + t_0(2)*30*24*3600 +t_0(1)*12*30*24*3600;
             [~, rho] = atmosnrlmsise00(h,lat,long,2017,153,UTC,'None');   %UM: [kg/m^3]
             rho = rho(6)*1e9;                                      %UM: [kg/km^3]

        elseif h > 11 && h < 86 %after troposphere
              h = h*1e3;
              [rho, ~, ~, ~, ~, ~] = atmos(h);      %UM: [kg/m^3]
              rho = rho*1e9;                        %UM: [kg/km^3]
        else 
              h = h*1e3;
              [rho, ~, ~, ~, ~]=tropos(h);          %UM: [kg/m^3]
              rho = rho*1e9;                        %UM: [kg/km^3]
        end
end