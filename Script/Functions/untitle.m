function [lat,long,hFig]=eci2ecf(varargin)
% La function è per definire la terna fissa con la terra, attraverso la
% rotazione della terna centrata inerziale (ECI)
%
% INPUTs  
% [x,y,z] 
%
% OUTPUTs
% lat : latitude
% long : longitude
% h fig : Earth figure handler
%
switch nargin
        case 2
            X0 = varargin{1}(1);
            Y0 = varargin{1}(2);
            Z0 = varargin{1}(3);
        case 1
            X0 = 0;
            Y0 = 0;
            Z0 = 0;
    otherwise
        error('Check the help of this function');
end
%Define Earth's topography surface
Earth = imread(['Textures\' earth '.jpg']);
% Define surface settings
        props.FaceColor= 'texture';
        props.EdgeColor = 'none';
        props.FaceLighting = 'gouraud';
        props.Cdata = Earth;



end
