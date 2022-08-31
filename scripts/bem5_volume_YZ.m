%   This script accurately computes and displays electric fields sampled on
%   a cross-section (sagittal plane) via the FMM method with accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Load/prepare data
planeABCD = [1 0 0 -X*1e-3];    % Equation of the plane of the cross-section (for neighbor triangle search speedup)
%%  Post processing parameters
component   = 4;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);

%%  Define dimensionless radius of the integration sphere (for precise integration near boundaries)
R = 4;

%%  Define observation points in the cross-section (MsxMs observation points)    
Ms = 500;
%   sagittal plane
y = linspace(ymin, ymax, Ms);
z = linspace(zmin, zmax, Ms);
[Y0, Z0]  = meshgrid(y, z);
clear pointsYZ;
pointsYZ(:, 1) = X*ones(1, Ms^2);
pointsYZ(:, 2) = reshape(Y0, 1, Ms^2);
pointsYZ(:, 3) = reshape(Z0, 1, Ms^2);

%   Set up enclosing tissues (optional)
pol = 1; % 1 - skin; 2 - skull; 3 - CSF
EofYZ_closed = close_meshpolygon(EofYZ{pol}, PofYZ{pol});
poly    = meshpolygon(PofYZ{pol}, EofYZ_closed);  % cross-section in the form of an oriented polygon
in      = inpolygon(pointsYZ(:, 2), pointsYZ(:, 3), poly(:, 2), poly(:, 3));

%  Assign tissue types to observation points (required for current density plot)
obsPointTissues = assign_tissue_type_volume(pointsYZ*1e-3, normals, Center, Indicator);

%% Find the E-feld at each observation point in the cross-section         
tic
pointsYZ       = 1e-3*pointsYZ;     % Convert back to m
Epri           = zeros(Ms*Ms, 3);
Esec           = zeros(Ms*Ms, 3);
Epri(in, :)    = bemf3_inc_field_electric(strcoil, pointsYZ(in, :), dIdt, mu0);      
Esec(in, :)    = bemf5_volume_field_electric(pointsYZ(in, :), c, P, t, Center, Area, normals, R, planeABCD);
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
figure;
%  E-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = abs(Etotal(:, component));
end
th1 = 120;          %   in V/m
th2 = 0;            %   in V/m
levels      = 20;
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('Distance y, mm');
ylabel('Distance z, mm');
title(['E-field (V/m), ', label, '-component in the sagittal plane.']);
 
% E-field plot: tissue boundaries
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% E-field plot: general settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([ymin ymax zmin zmax]);
grid on; set(gcf,'Color','White');


%% Plot the current distribution in the cross-section
figure;
%  J-field plot: contour plot
if component == 4
    temp      = abs(sqrt(dot(Jtotal, Jtotal, 2)));
else
    temp      = abs(Itotal(:, component));
end
th1 = 120;          %   in V/m
th2 = 0;            %   in V/m
levels      = 20;
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('Distance y, mm');
ylabel('Distance z, mm');
title(['Current (A/m^2), ', label, '-component in the sagittal plane.']);

% J-field plot: tissue boundaries
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
% J-field plot: General settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([ymin ymax zmin zmax]);
grid on; set(gcf,'Color','White');

%% Verification: Plot obs pt tissue assignments against the contours
% The tissue assignment may fail for bad meshes or complicated geometry.
% This plot will show which tissue types were used for the prior plot
color = prism(length(tissue)); color(4,:) = [0 1 1];
figure;
legendEntries = {};
for j = 0:length(tissue)
    if(j == 0)
        c_temp = [0 0 0];
    else
        c_temp = color(j,:);
    end
    ptsY = pointsYZ(obsPointTissues == j, 2) * 1e3;
    ptsZ = pointsYZ(obsPointTissues == j, 3) * 1e3;
    scatter(ptsY, ptsZ, 4, c_temp, 'filled'); hold on;
end
% Display the tissue boundaries
color = ones(length(tissue), 3);
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
%   General settings 
axis 'equal';  axis 'tight';     
axis([ymin ymax zmin zmax]);
grid on; set(gcf,'Color','White');
title('Tissue types assigned to observation points');
xlabel('Distance y, mm');
ylabel('Distance z, mm');
