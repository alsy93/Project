function [vargout] = flightEnvelope(h,v,par)

% This function provides the flight envelope constraints

% 1)Equilibrium trajectory [gamma is zero and also the dgamma/dt is zero]:
%   maximum lift ceiling

    h_eq = log((1.225*1e9.*v.^2)./(2*par.alpha.*(par.g - v.^2./par.Re)))./...
           (900/par.Re - 1/par.Re .* (2*par.g - (v.^2./par.Re)./(par.g - (v.^2./par.Re)))); 
% 2)Dynamic pressure limit


%     function f_eq = findheq(v,par)
%      
%         
%         f_eq =@(h_eq)  h_eq - 1./(par.g./(v.^2) - rho./(2*par.alpha)) + par.Re;
%         rho = varrho(h_eq);
%     end

    %v = v*;
    %g = g(h)
%     
%     nV = length(v);
%     h_eq = zeros(nV,1);
%     f_eq =@(h_eq,v)  h_eq - 1./(par.g./(v.^2) - varrho(h_eq)./(2*par.alpha)) + par.Re; 
%     myOptions = optimoptions('fsolve','Display','off');
%     for i=1:nV
%         fprintf(['Run: ',num2str(i),'of ',num2str(nV),'\n']);
%         tempFun=@(xx) f_eq(xx,v(i));
%         [h_eq(i),~] = fsolve(tempFun,100,myOptions);
%     end
%     
%     clear tempFun i 
%     options = optimoptions('fsolve','OutputFcn',@myOutFcn,'Diagnostics','on','Display','iter-detailed');


% Plot solutions




    figure();hold on; grid on; grid minor
    title('System and operational constraints');
    xlabel('Relative velocity');
    ylabel('Altitude');
    legend();
    vargout(1) = plot (v,h,'-b','LineWidth',2,'DisplayName','Computed solution');
    vargout(2) = plot(v,h_eq,'-r','LineWidth',2,'DisplayName','Maximum lift ceiling');
    
end



    
    

