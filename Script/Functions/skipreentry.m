function [ts,ys] = skipreentry(y0,time,par)
        
    [t1, y1] = integratorskip(y0,time,par);

    N = length(y1);
    y0_skip = [y1(N,1), y1(N,2), y1(N,3), y1(N,4), y1(N,5), y1(N,6)];
    t0_skip = [t1(N) 10000];
    
    [t, y] = integratorskip1(y0_skip,t0_skip,par);
    
    N1 = length(y);
    y1_skip = [y(N1,1), y(N1,2), y(N1,3), y(N1,4), y(N1,5), y(N1,6)];
    t1_skip = [t(N1) 10000];
    
    [t2, y2] = integratorskip2(y1_skip,t1_skip,par);
   
    
    ts = [t1;t; t2];
    ys =[y1;y; y2];
    
    
    %Plotting results
     
    figure()
    
        ax1 = subplot(2,3,1);
        plot(ax1,ts,ys(:,1))
        title('Velocity')
        grid on;grid minor
        
        ax2 = subplot(2,3,2);
        plot(ax2,ts,ys(:,2))
        title('Flight path angle')
        grid on;grid minor
        
        ax3 = subplot(2,3,3);
        plot(ax3,ts,ys(:,3))
        title('Altitude')
        grid on;grid minor
        
        ax4 = subplot(2,3,4);
        plot(ax4,ts,ys(:,4))
        title('Latitude')
        grid on;grid minor
        
        ax5 = subplot(2,3,5);
        plot(ax5,ts,ys(:,5))
        title('Longitude')
        grid on;grid minor
        
        ax6 = subplot(2,3,6);
        plot(ax6,ts,ys(:,6))
        title('Heading angle')
        grid on;grid minor
      
  
        

end