function [vargout] = groundtrack(lat,long)

%Graphic init
    handler_fig = figure();
    handler_fig.Name = 'Ground Track Plot';
    grid on
    hold on
%Load the image of earth surface
    immagine = imread('Textures\earth.jpg');
%Resize in a proper way the graphic
    imagesc([-180 180],[-90 90],flip(immagine))
    axis([0 100 0 80]);
    
    vargout(1) = plot(long,lat,'Or');

end

