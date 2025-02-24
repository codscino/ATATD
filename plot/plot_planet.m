function h = plot_planet(x0, y0, z0, radius, planet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_planet - Plot a textured planet sphere in 3D space
%
%   Inputs:
%       x0, y0, z0 - Coordinates of the planet's center.
%       radius     - Radius of the planet.
%       planet     - Filename of the image file to be used as the
%                    planet's texture.
%
%   Outputs:
%       h          - Handle to the plotted surface object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Read the texture image for the planet
    planetTexture = imread(planet); 
    planetTexture = flipud(planetTexture);

    % Create a sphere (by default 21-by-21 grid) and scale it by the radius
    [x, y, z] = sphere;
    baseX = radius * x;  
    baseY = radius * y;  
    baseZ = radius * z;  
    
    % Shift the sphere to the desired center (x0, y0, z0)
    x = baseX + x0;
    y = baseY + y0;
    z = baseZ + z0;
    
    % Plot the surface with the texture map and return its handle
    h = surf(x, y, z, 'FaceColor', 'texturemap', ...
                  'CData', planetTexture, 'EdgeColor', 'none');
              
    % Save the base (unshifted) coordinates and radius in the object's UserData
    h.UserData.baseX = baseX;
    h.UserData.baseY = baseY;
    h.UserData.baseZ = baseZ;
    h.UserData.radius = radius;
end