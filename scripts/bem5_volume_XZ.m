%   This script accurately computes and displays electric fields sampled on
%   a cross-section (coronal plane) via the FMM method with accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Load/prepare data
planeABCD = [0 1 0 -Y*1e-3];    % Cross section plane (for neighbor triangle search acceleration)
%%  Post processing parameters
component   = 4;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);

%%  Define dimensionless radius of the integration sphere (for precise integration near boundaries)
R = 4;

%%  Define observation points in the cross-section (MsxMs nodal points)    
Ms = 200;
%   Coronal plane
x = linspace(xmin, xmax, Ms);
z = linspace(zmin, zmax, Ms);
[X0, Z0]  = meshgrid(x, z);
clear pointsXZ;
pointsXZ(:, 1) = reshape(X0, 1, Ms^2);
pointsXZ(:, 2) = Y*ones(1, Ms^2);
pointsXZ(:, 3) = reshape(Z0, 1, Ms^2);  

%   Set up enclosing tissues (optional)
pol = 1; % 1 - skin; 2 - skull; 3 - CSF
EofXZ_closed = close_meshpolygon(EofXZ{pol}, PofXZ{pol});
poly    = meshpolygon(PofXZ{pol}, EofXZ_closed);  % cross-section in the form of an oriented polygon
in      = inpolygon(pointsXZ(:, 1), pointsXZ(:, 3), poly(:, 1), poly(:, 3));

%  Assign tissue types to observation points (required for current density plot)
obsPointTissues = assign_tissue_type_volume(pointsXZ*1e-3, normals, Center, Indicator);

%%  Find the E-field at each observation point in the cross-section (MsxMs nodal points)
tic
pointsXZ       = 1e-3*pointsXZ;     % Convert back to m
Epri           = zeros(Ms*Ms, 3);
Esec           = zeros(Ms*Ms, 3);
Epri(in, :)    = bemf3_inc_field_electric(strcoil, pointsXZ(in, :), dIdt, mu0);      
Esec(in, :)    = bemf5_volume_field_electric(pointsXZ(in, :), c, P, t, Center, Area, normals, R, planeABCD);
Etotal         = Epri + Esec;   
fieldPlaneTime = toc  

%% Calculate current density at each observation point
% Observation points in free space were originally assigned tissue code 0, 
% so provide an extra "Free Space" conductivity entry and point free space obs pts to that entry.
condTemp = [cond 0];                                        
obsPointTissues(obsPointTissues == 0) = length(condTemp);

% Calculate current density J = sigma * E
condTempExpanded = transpose(condTemp(obsPointTissues));
Jtotal = repmat(condTempExpanded, 1, 3).*Etotal;

%%  Plot the E-field in the cross-section
%  E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = abs(Etotal(:, component));
end
th1 = 120;          %   in V/m
th2 = 0;            %   in V/m
levels      = 20;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, z);
xlabel('Distance x, mm');
ylabel('Distance z, mm');
title(['E-field (V/m), ', label, '-component in the coronal plane.']);
 
% E-field plot: tissue boundaries
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([xmin xmax zmin zmax]);
grid on; set(gcf,'Color','White');


%% Plot the current distribution in the cross-section
figure;
% J-field plot: Contour plot
if component == 4
    temp      = abs(sqrt(dot(Jtotal, Jtotal, 2)));
else
    temp      = abs(Jtotal(:, component));
end
th1 = 120;          %   in V/m
th2 = 0;            %   in V/m
levels      = 20;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, z);
xlabel('Distance x, mm');
ylabel('Distance z, mm');
title(['Current density (A/m^2), ', label, '-component in the coronal plane.']);

% J-field plot: tissue boundaries
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXY
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% J-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([xmin xmax zmin zmax]);
grid on; set(gcf,'Color','White');


%% Verification: Plot obs pt tissue assignments against the contours
% The tissue assignment may fail for bad meshes or complicated geometry.
% This plot will show which tissue types were used for the prior plot
color = prism(length(tissue)); color(4,:) = [0 1 1];
figure;
for j = 0:length(tissue)
    if(j == 0)
        c_temp = [0 0 0];
    else
        c_temp = color(j,:);
    end
    ptsX = pointsXZ(obsPointTissues == j, 1) * 1e3;
    ptsZ = pointsXZ(obsPointTissues == j, 3) * 1e3;
    scatter(ptsX, ptsZ, 4, c_temp, 'filled'); hold on;
end
% Display the tissue boundaries
color = ones(length(tissue), 3);
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
%   General settings 
axis 'equal';  axis 'tight';     
axis([xmin xmax zmin zmax]);
grid on; set(gcf,'Color','White');
title('Tissue types assigned to observation points');
xlabel('Distance x, mm');
ylabel('Distance z, mm');