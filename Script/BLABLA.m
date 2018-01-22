function [tt,yy,collt,collX] = BLABLA()

collX = [];
collt = [];
 


myOptions = odeset('OutputFcn',@myOutFcn);
[tt,yy] = ode45(@myFun,[0,10],5,myOptions);


function dx = myFun(t,x)
    dx = -2*x;
end

function status = myOutFcn(t,x,flag)
    switch flag
        case 'init'
            collX = [collX,x(1)];
            collt = [collt,t(1)];
        case []
            collX = [collX,x];
            collt = [collt,t];
        otherwise
            
    end
    
    status = 0;
end

end


