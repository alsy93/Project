function [vargout] = GroundTrack(v_ECI,lat,long)

%Graphic init
    handler_fig = figure();
    handler_fig.Name = 'Ground Track Plot';
    handler_fig.Position = [screenSize(3)/4 screenSize(4)/4 screenSize(3)/2 screenSize(4)/2];
    grid on
    hold on
%Load the image of earth surface
    immagine = imread('Textures\earth.jpg');
%Resize in a proper way the graphic
    imagesc([-180 180],[-90 90],flip(immagine))
    axis([-180 180 -90 90]);
    
   plot(lat,long);

end

