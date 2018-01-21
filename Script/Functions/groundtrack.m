function [vargout] = groundtrack(long,lat)

%Graphic init
    handler_fig = figure();
    handler_fig.Name = 'Ground Track Plot';
    grid on
    hold on
%Load the image of earth surface
    immagine = imread('Textures\earth.jpg');
%Resize in a proper way the graphic
    imagesc([-180 180],[-90 90],flip(immagine))
    axis([0 90 0 90]);
    legend();
     
    
    vargout(1) = plot(long,lat,'Or','DisplayName','Reentry path');
    vargout(2) = plot(long(1),lat(1),'Og','DisplayName','Departure');
    n = length(long);
    vargout(3) = plot(long(n),lat(n),'Xy','LineWidth',2,'DisplayName','Arrival');
    

end

