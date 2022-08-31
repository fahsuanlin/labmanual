%   This script computes and plots the incident magnetic field magnitude
%   (or any of the components) for any brain compartment surface/interface
%   (plots the surface field + optionally coil geometry) 
%
%   Copyright SNM/WAW 2017-2020

%%   Compute the B-field for all surfaces/interfaces
tissue_to_plot = 'GM';

objectnumber = find(strcmp(tissue, tissue_to_plot));              
Bmag         = sqrt(dot(Binc(Indicator==objectnumber, :), Binc(Indicator==objectnumber, :), 2));

%%   Graphics
figure;
bemf2_graphics_surf_field(P, t, Bmag, Indicator, objectnumber);
title(strcat('Solution: B-field magnitude in T for: ', tissue{objectnumber}));

% Coil centerline graphics    
bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
hold on;
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 3);

% General
axis tight; view(30, 35); camzoom(1);
brighten(0.4);