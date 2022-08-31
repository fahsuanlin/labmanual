%   This script plots the electric potential of the secondary field for
%   any brain compartment surface (plots the density + optionally
%   coil geometry).
%
%   Copyright SNM/WAW 2017-2020

%%   Graphics
tissue_to_plot = 'GM';

objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = Pot(Indicator==objectnumber);
figure;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Electric potential of the secondary field for: ', tissue{objectnumber}));

% Coil centerline graphics 
bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
hold on;
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

% General
axis tight; view(30, 35), camzoom(1);

brighten(0.4);