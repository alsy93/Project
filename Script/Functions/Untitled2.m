h = [35; 36; 37; 38];
g = [60;61;62; 63];
v0 = [5; 6; 7; 8];
for i = 1:4
    v_body(:,i) = [v0(i)*cosd(h(i)) + v0(i)*cosd(g(i));...
                   v0(i)*sind(h(i)); v0(i)*sind(g(i))];
R1=[cosd(h(i)) -sind(h(i)) 0; sind(h(i)) cosd(h(i)) 0;0 0 1];
v(:,i) = R1*v_body(:,i); 
end

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


