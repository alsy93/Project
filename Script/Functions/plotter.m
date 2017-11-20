function plothandler = plotter(Name,varargin)

%Questa function serve per disegnare la terra o il centro di pressione o 
%baricentro della Soyuz.

% INPUTs
% Name : 'earth', 'Cg', 'Cp'  
% [X0,Y0,Z0]: position of the object (for earth use only the name as the
%             input)
%

% OUTPUTs
% planet_handler      planet or object figure handler
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

notFound = 0;
switch Name
    case {'earth',1}
            Name = 'earth';
            Radius=0.63781600000000E+04; % [km] From DITAN
    case {'Cg',2}
            Name='Cg';
            Radius=1; %[km]
    case {'Cp',3}
            Name='Cp';
            Radius=1; %[km]
    otherwise
            disp('This function can''t rapresent this object: check the list of supported one')
            notFound = 1;
end

if notFound == 0
   
    %Load and define topographic data

    object = imread(['Textures\' Name '.jpg']);

    % Define surface settings

            props.FaceColor= 'texture';
            props.EdgeColor = 'none';
            props.FaceLighting = 'gouraud';
            props.Cdata = object;

    % Create the sphere with Earth or object topography and adjust colormap

        
        [X, Y, Z] = ellipsoid(X0,Y0,Z0,Radius,Radius,Radius,200);
        plothandler = surface(X,Y,flip(Z),props);
end
return
